#### selection_rate_analysis.R by Jack Common
#### Code to generate figures and models that test the hypotheses that phage diversity
#### will a) negate the benefits to the host of CRISPR diversity and b) that host diversity
#### will dilute individual susceptible genotypes.

rm(list=ls())
options(show.signif.stars = F)

#### ---- Dependencies ---- ####
library(scales)
library(ggstatsplot)
library(cowplot)
library(ggplot2)
library(plyr)
library(lme4)
library(magrittr)
library(tidyverse)
library(lmerTest)
library(MCMCglmm)

#### ---- Functions ---- ####
## Compare AIC values for model fit
# This function extracts the AIC for each GLM, and then compares the absolute relative 
# differences for each AIC to the model with the lowest AIC. This acts as a measure of 
# model fit. More can be found at: http://faculty.washington.edu/skalski/classes/QERM597/papers_xtra/Burnham%20and%20Anderson.pdf

compare_AICs = function(df){          # df is a dataframe of AIC values 
  print(df)                           # prints the origina AIC values 
  col_len = length(df[,2])            # extracts the number of number of models
  AIC_min = abs(min(df[,2]))          # finds the minimum AIC value
  for (i in seq(1, col_len, 1)){      # loop through the AIC values and prints the absolute differences from AIC_min
    print( (abs(df[i,2])) - AIC_min)
  }
}

## Convert log-odds to probabilities from binomial GLM output
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}
treatment_names <- list(
  "3" = "3-clone",
  "6" = "6-clone",
  "12" = "12-clone",
  "24" = "24-clone"
)

treatment_labeller <- function(variable, value) {
  return(treatment_names[value])
}

timepoint_names <- list(
  "1" = "1 dpi",
  "2" = "2 dpi",
  "3" = "3 dpi"
)

timepoint_labeller <- function(variable, value) {
  return(timepoint_names[value])
}

treatment_list <- c("3-clone", "6-clone", "12-clone", "24-clone")

pfu_labels <- c(
  expression("0"), 
  expression('10'^2*''), 
  expression('10'^4*''),
  expression('10'^6*''),
  expression('10'^8*''),
  expression('10'^10*''),
  expression("10"^12*"")
)


#### ---- Original data ---- ####
rates <- read.csv("./Experiment 1/original_data/reldiv_exp1_master.csv", header=T, skip = 1)  %>% 
  select(Treatment, Timepoint, Replicate, BIM_mix, Escape_phage, Tracked_BIM,
         rCRISPR, rBIM) %>% 
  gather("rCRISPR", "rBIM", key="Strain", value="r", factor_key = T) %>% 
  filter(Timepoint %in% c("1", "2", "3"), r!='#NUM!')

rates$Timepoint %<>% as.factor()
rates$Replicate %<>% as.factor
rates$Tracked_BIM %<>% as.factor()
rates$Escape_phage %<>% as.factor
rates$Treatment %<>% as.factor
rates$r %<>% as.numeric()

# Relevel the datasets so things make sense
rates$Treatment %<>% relevel(ref="24")
rates$Treatment %<>% relevel(ref="12")
rates$Treatment %<>% relevel(ref="6")
rates$Treatment %<>% relevel(ref="3")

#### ---- CRISPR raw data boxplot ---- ####

CRISPR_rate_plot <- ggplot(aes(x=Treatment, y=r), 
                               data=filter(rates, Strain=="rCRISPR", Timepoint!="0"))+
  geom_boxplot(width=0.5, outlier.shape = NA)+
  geom_jitter(width=0.2, colour="black", alpha=0.5, pch=21)+
  geom_hline(yintercept=0, linetype=2, size=1)+
  #stat_summary(fun.y="mean", geom="point", colour="red")+
  facet_wrap(~Timepoint, labeller = timepoint_labeller, ncol=2)+
  labs(y="Selection rate (day -1) of CRISPR", x="CRISPR x phage diversity")+
  #scale_x_discrete(labels = treatment_list)+
  theme_cowplot()+
  theme(axis.title = element_text(face="bold", size=16),
        panel.grid.major = element_line(colour="lightgrey", size=.25),
        strip.text = element_text(face="bold"),
        axis.text.x = element_text(size=10))+
  #coord_cartesian(ylim=c(-2.5, 2.5))+ #use this to "zoom in" on values closer to zero
  NULL

last_plot()

ggsave("rate_CRISPR_raw.png", CRISPR_rate_plot, path="./Experiment 1/docs/", device="png", dpi=600,
       width=28, height = 20, units = c("cm"))

#### ---- BIM raw data boxplot ---- ####
BIM_rate_plot <- ggplot(aes(x=Treatment, y=r), 
                            data=filter(rates, Strain=="rBIM", Timepoint!="0"))+
  geom_boxplot(width=0.5, outlier.shape = NA)+
  geom_jitter(width=0.2, colour="black", alpha=0.5, pch=21)+
  geom_hline(yintercept=0, linetype=2, size=1)+
  #stat_summary(fun.y="mean", geom="point", colour="red")+
  labs(y="Selection rate (day -1) of\nsusceptible clone", x="CRISPR x phage diversity")+
  #scale_x_discrete(labels= treatment_list[3:7])+
  theme_cowplot()+
  #facet_wrap(~Timepoint, scales="free")+
  facet_wrap(~Timepoint, labeller=timepoint_labeller, ncol=2)+
  theme(axis.title = element_text(face="bold", size=16),
        panel.grid.major = element_line(colour="lightgrey", size=.25),
        strip.text = element_text(face="bold"),
        axis.text.x = element_text(size=10))+
  NULL
last_plot()

ggsave("rate_BIM_raw.png", BIM_rate_plot, path="./Experiment 1/docs/", device="png", dpi=600,
       width=28, height = 20, units = c("cm"))

### ---- lme4 CRISPR models ---- ####
d.CRISPR <- rates %>% 
  na.exclude %>% 
  filter(Strain=="rCRISPR", Timepoint!="0")

d.CRISPR$Treatment %<>% relevel(ref="24")
d.CRISPR$Treatment %<>% relevel(ref="12")
d.CRISPR$Treatment %<>% relevel(ref="6")
d.CRISPR$Treatment %<>% relevel(ref="3")

d.CRISPR$Timepoint %<>% relevel(ref="3")
d.CRISPR$Timepoint %<>% relevel(ref="2")
d.CRISPR$Timepoint %<>% relevel(ref="1")

m1 <- lmer(r~Treatment + (1|Replicate), data=d.CRISPR)
m2 <- lmer(r~Timepoint + (1|Replicate), data=d.CRISPR)
m3 <- lmer(r~Treatment + Timepoint +(1|Replicate), data=d.CRISPR)
m4 <- lmer(r~Treatment * Timepoint + (1|Replicate), data=d.CRISPR)

AIC(m1, m2, m3, m4) %>% 
  compare_AICs()

plot(m1)
plot(m2)
plot(m3)
plot(m4)

summary(m3)
anova(m3, type="marginal")
ggcoefstats(m2)

predict(m3, d.CRISPR, type="response") %>% 
  unique

### ---- lme4 BIM models ---- ####
d.BIM <- rates %>% 
  na.exclude %>% 
  filter(Strain=="rBIM", Timepoint!="0")

d.BIM$Treatment %<>% relevel(ref="24")
d.BIM$Treatment %<>% relevel(ref="12")
d.BIM$Treatment %<>% relevel(ref="6")
d.BIM$Treatment %<>% relevel(ref="3")

d.BIM$Timepoint %<>% relevel(ref="3")
d.BIM$Timepoint %<>% relevel(ref="2")
d.BIM$Timepoint %<>% relevel(ref="1")

m1 <- lmer(r~Treatment + (1|Replicate), data=d.BIM)
m2 <- lmer(r~Timepoint + (1|Replicate), data=d.BIM)
m3 <- lmer(r~Treatment + Timepoint +(1|Replicate), data=d.BIM)
m4 <- lmer(r~Treatment * Timepoint + (1|Replicate), data=d.BIM)

AIC(m1, m2, m3, m4) %>% 
  compare_AICs()

plot(m1)
plot(m2)
plot(m3)
plot(m4)

summary(m4)
summary(m3)
anova(m4, type="marginal")
ggcoefstats(m2)

predict(m2, d.CRISPR, type="response") %>% 
  unique


