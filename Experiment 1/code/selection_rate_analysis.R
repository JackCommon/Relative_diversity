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

#### ---- CRISPR boxplot ---- ####

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

#### ---- BIM boxplot ---- ####
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


#### ---- MCMC CRISPR models ---- ####
library(data.table)
library(ape)
library(coda)
library(lattice)
library(VCVglmm)
library(MCMCglmm)
library(aod)

d.CRISPR <- d %>% 
  na.exclude %>% 
  filter(Strain=="rCRISPR")

d.CRISPR$Treatment %<>% relevel(ref="24-clone_control")
d.CRISPR$Treatment %<>% relevel(ref="1-clone_control")
d.CRISPR$Treatment %<>% relevel(ref="24-clone")
d.CRISPR$Treatment %<>% relevel(ref="12-clone")
d.CRISPR$Treatment %<>% relevel(ref="6-clone")
d.CRISPR$Treatment %<>% relevel(ref="3-clone")
d.CRISPR$Treatment %<>% relevel(ref="1-clone")

d.CRISPR$Timepoint %<>% relevel(ref="3")
d.CRISPR$Timepoint %<>% relevel(ref="2")
d.CRISPR$Timepoint %<>% relevel(ref="1")
d.CRISPR$Timepoint %<>% relevel(ref="0")

m1 <- MCMCglmm(r~Treatment, random=~Replicate, 
               data=d.CRISPR,
               nitt = 13000, burnin = 3000, thin = 10, pr = TRUE)
m2 <- MCMCglmm(r~Timepoint, random=~Replicate, 
               family="gaussian", data=filter(d.CRISPR, Escape_phage!="ancestral"),
               nitt = 13000, burnin = 3000, thin = 10, pr = TRUE)
m3 <- MCMCglmm(r~Treatment + Timepoint , random=~Replicate, 
               family="gaussian", data=filter(d.CRISPR,Treatment!="1-clone"),
               nitt = 13000, burnin = 3000, thin = 10, pr = TRUE)
m4 <- MCMCglmm(r~Treatment * Timepoint, random=~Replicate,
               family="gaussian", data=d.CRISPR,
               nitt = 13000, burnin = 3000, thin = 10, pr = TRUE)

m1$DIC; m2$DIC; m3$DIC; m4$DIC

plot(m1)
plot(m2)
plot(m3)

Wald.test.auto(m1)
Wald.test.auto(m2)
Wald.test.auto(m3)

summary(m3)

MCMCfixplot(m3)

#### Coefficient plot: CRISPR ####
# Set up target dataframe to store modes and HPDs
m3 <- MCMCglmm(r~Treatment + Timepoint , random=~Replicate, 
               family="gaussian", data=d.CRISPR,
               nitt = 13000, burnin = 3000, thin = 10, pr = TRUE)

CRISPR.coefs <- data.frame(term = factor(9), 
                           mode = numeric(9), 
                           l.95 = numeric(9), h.95 = numeric(9), 
                           l.89 = numeric(9), h.89 = numeric(9), 
                           l.67 = numeric(9), h.67 = numeric(9))

CRISPR.coefs$term <- c("1-clone", "3-clone", "6-clone", "12-clone",
                       "24-clone", '1-clone\n(ancestral phage)', "24-clone\n(ancestral phage)",
                       "2 dpi", "3 dpi") %>% 
  as.factor

# Get posterior coefficients
CRISPR.coefs$mode <- posterior.mode(m3$Sol)[1:9]
CRISPR.coefs$l.95 <- HPDinterval(m3$Sol, prob=0.95) %>% 
  as.data.frame() %>%    
  slice(1:9) %>% 
  select(lower) %>% 
  .[1:9,1]
CRISPR.coefs$h.95 <- HPDinterval(m3$Sol, prob=0.95) %>% 
  as.data.frame() %>%    
  slice(1:9) %>% 
  select(upper) %>% 
  .[1:9,1]

CRISPR.coefs$l.89 <- HPDinterval(m3$Sol, prob=0.89) %>% 
  as.data.frame() %>%    
  slice(1:9) %>% 
  select(lower ) %>% 
  .[1:9,1]
CRISPR.coefs$h.89 <- HPDinterval(m3$Sol, prob=0.89) %>% 
  as.data.frame() %>%    
  slice(1:9) %>% 
  select(upper) %>% 
  .[1:9,1]

CRISPR.coefs$l.67 <- HPDinterval(m3$Sol, prob=0.67) %>% 
  as.data.frame() %>%    
  slice(1:9) %>% 
  select(lower) %>% 
  .[1:9,1]
CRISPR.coefs$h.67 <- HPDinterval(m3$Sol, prob=0.67) %>% 
  as.data.frame() %>%    
  slice(1:9) %>% 
  select(upper) %>% 
  .[1:9,1]

CRISPR.coefs$term %<>% relevel(ref="3 dpi")
CRISPR.coefs$term %<>% relevel(ref="2 dpi")
CRISPR.coefs$term %<>% relevel(ref="24-clone\n(ancestral phage)")
CRISPR.coefs$term %<>% relevel(ref="24-clone")
CRISPR.coefs$term %<>% relevel(ref="12-clone")
CRISPR.coefs$term %<>% relevel(ref="6-clone")
CRISPR.coefs$term %<>% relevel(ref="3-clone")
CRISPR.coefs$term %<>% relevel(ref="1-clone\n(ancestral phage)")
CRISPR.coefs$term %<>% relevel(ref="1-clone")

CRISPR.coefs <- CRISPR.coefs %>% 
  filter(term!="2 dpi", term!="3 dpi")


pd1 <- position_dodge(1)
pd2 <- position_dodge(0.7)
p1 <- ggplot(aes(y=mode, x=term), data=CRISPR.coefs)+
  geom_errorbar(aes(ymin=l.67, ymax=h.67), width=0, size=4, alpha=0.5)+
  geom_errorbar(aes(ymin=l.89, ymax=h.89), width=0, size=2, alpha=0.5)+
  geom_errorbar(aes(ymin=l.95, ymax=h.95), width=0)+
  geom_point(fill="white", pch=21, colour="black", size=2)+
  geom_hline(yintercept=0, linetype=2)+
  #coord_flip()+
  cowplot::theme_cowplot()+
  coord_cartesian(ylim=c(-2.592837,  2.748592))+
  labs(y="CRISPR selection rate", x="CRISPR allele diversity")+ 
  #scale_y_continuous(breaks=seq(-0.5, 1.5, 0.5))+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(face="bold", size=16))+
  NULL
p1

#

### Plotting predictions: CRISPR ####

d.CRISPR$m3.preds <- predict(m3, d.CRISPR, type="response")
d.CRISPR$m4.means <- predict(m4, d.CRISPR, type="response")

d.CRISPR$m4.95.lower <- predict(m4, d.CRISPR, type="response",
                                interval="confidence", level=0.95)[,"lwr"]
d.CRISPR$m4.95.upper <- predict(m4, d.CRISPR, type="response",
                                interval="confidence", level=0.95)[,"upr"]
d.CRISPR$m4.89.lower <- predict(m4, d.CRISPR, type="response",
                                interval="confidence", level=0.89)[,"lwr"]
d.CRISPR$m4.89.upper <- predict(m4, d.CRISPR, type="response",
                                interval="confidence", level=0.89)[,"upr"]
d.CRISPR$m4.67.lower <- predict(m4, d.CRISPR, type="response",
                                interval="confidence", level=0.67)[,"lwr"]
d.CRISPR$m4.67.upper <- predict(m4, d.CRISPR, type="response",
                                interval="confidence", level=0.67)[,"upr"]

CRISPR.preds <- d.CRISPR %>% 
  select(Treatment, Timepoint, m4.means, 
         m4.95.lower, m4.95.upper,
         m4.89.lower, m4.89.upper,
         m4.67.lower, m4.67.upper) %>% 
  unique %>% 
  filter(Timepoint!="0")

pd <- position_dodge(0.5)

p2 <- ggplot(aes(x=Treatment, y=m4.means), data=CRISPR.preds)+
  labs(x="Days post-infection", y="CRISPR selection rate")+
  geom_hline(yintercept = 0, linetype=2, colour="black")+
  geom_errorbar(aes(ymin=m4.67.lower, ymax=m4.67.upper),
                width=0, size=4, alpha=0.5, position=pd)+
  geom_errorbar(aes(ymin=m4.89.lower, ymax=m4.89.upper),
                width=0, size=2, alpha=0.5, position=pd)+
  geom_errorbar(aes(ymin=m4.95.lower, ymax=m4.95.upper), 
                width=0, position=pd)+
  geom_point(aes(position=Treatment), pch=21, fill="white", 
             colour="black", position=pd, size=2)+
  facet_wrap(~Timepoint)+
  cowplot::theme_cowplot()+
  scale_y_continuous(breaks=seq(-4,2,1))+
  coord_cartesian(ylim=c(-4.226857, 2.425297))+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(face="bold", size=16),
        #panel.grid.major = element_line(colour="lightgrey", size=.25),
        strip.text = element_text(face="bold"))+
  NULL

p2


#### ---- BIM models ---- ####

# d.BIM <- d %>% 
#   na.exclude %>% 
#   filter(Strain=="rBIM", Timepoint!="0",
#          Treatment%in%c("3-clone", "6-clone", "12-clone",
#                         "24-clone", "1-clone_control", "24-clone_control"))

d.BIM <- d %>% 
  na.exclude %>% 
  filter(Strain=="rBIM", Timepoint!="0")

d.BIM$Treatment %<>% relevel(ref="24-clone_control")
d.BIM$Treatment %<>% relevel(ref="1-clone_control")
d.BIM$Treatment %<>% relevel(ref="24-clone")
d.BIM$Treatment %<>% relevel(ref="12-clone")
d.BIM$Treatment %<>% relevel(ref="6-clone")
d.BIM$Treatment %<>% relevel(ref="3-clone")
d.BIM$Treatment %<>% relevel(ref="1-clone")

d.BIM$Timepoint %<>% relevel(ref="3")
d.BIM$Timepoint %<>% relevel(ref="2")
d.BIM$Timepoint %<>% relevel(ref="1")


m1 <- MCMCglmm(r~Treatment, random=~Replicate, 
               data=d.BIM,
               nitt = 13000, burnin = 3000, thin = 10, pr = TRUE)
m2 <- MCMCglmm(r~Timepoint, random=~Replicate, 
               family="gaussian", data=d.BIM,
               nitt = 13000, burnin = 3000, thin = 10, pr = TRUE)
m3 <- MCMCglmm(r~Treatment + Timepoint , random=~Replicate, 
               family="gaussian", data=d.BIM,
               nitt = 13000, burnin = 3000, thin = 10, pr = TRUE)
m4 <- MCMCglmm(r~Treatment * Timepoint, random=~Replicate,
               family="gaussian", data=d.BIM,
               nitt = 13000, burnin = 3000, thin = 10, pr = TRUE)

m1$DIC; m2$DIC; m3$DIC; m4$DIC

plot(m1)
plot(m2)
plot(m3)

Wald.test.auto(m1)
Wald.test.auto(m2)
Wald.test.auto(m3)

summary(m3)

MCMCfixplot(m3)

#### Coefficient plot: BIMs ####
# Set up target dataframe to store modes and HPDs
m3 <- MCMCglmm(r~Treatment + Timepoint , random=~Replicate, 
               family="gaussian", data=d.BIM,
               nitt = 13000, burnin = 3000, thin = 10, pr = TRUE)

BIM.coefs <- data.frame(term = factor(9), 
                        mode = numeric(9), 
                        l.95 = numeric(9), h.95 = numeric(9), 
                        l.89 = numeric(9), h.89 = numeric(9), 
                        l.67 = numeric(9), h.67 = numeric(9))

BIM.coefs$term <- c("1-clone", "3-clone", "6-clone", "12-clone",
                    "24-clone", '1-clone\n(ancestral phage)', "24-clone\n(ancestral phage)",
                    "2 dpi", "3 dpi") %>% 
  as.factor

# Get posterior coefficients
BIM.coefs$mode <- posterior.mode(m3$Sol)[1:9]
BIM.coefs$l.95 <- HPDinterval(m3$Sol, prob=0.95) %>% 
  as.data.frame() %>%    
  slice(1:9) %>% 
  select(lower) %>% 
  .[1:9,1]
BIM.coefs$h.95 <- HPDinterval(m3$Sol, prob=0.95) %>% 
  as.data.frame() %>%    
  slice(1:9) %>% 
  select(upper) %>% 
  .[1:9,1]

BIM.coefs$l.89 <- HPDinterval(m3$Sol, prob=0.89) %>% 
  as.data.frame() %>%    
  slice(1:9) %>% 
  select(lower ) %>% 
  .[1:9,1]
BIM.coefs$h.89 <- HPDinterval(m3$Sol, prob=0.89) %>% 
  as.data.frame() %>%    
  slice(1:9) %>% 
  select(upper) %>% 
  .[1:9,1]

BIM.coefs$l.67 <- HPDinterval(m3$Sol, prob=0.67) %>% 
  as.data.frame() %>%    
  slice(1:9) %>% 
  select(lower) %>% 
  .[1:9,1]
BIM.coefs$h.67 <- HPDinterval(m3$Sol, prob=0.67) %>% 
  as.data.frame() %>%    
  slice(1:9) %>% 
  select(upper) %>% 
  .[1:9,1]

BIM.coefs$term %<>% relevel(ref="3 dpi")
BIM.coefs$term %<>% relevel(ref="2 dpi")
BIM.coefs$term %<>% relevel(ref="24-clone\n(ancestral phage)")
BIM.coefs$term %<>% relevel(ref="24-clone")
BIM.coefs$term %<>% relevel(ref="12-clone")
BIM.coefs$term %<>% relevel(ref="6-clone")
BIM.coefs$term %<>% relevel(ref="3-clone")
BIM.coefs$term %<>% relevel(ref="1-clone\n(ancestral phage)")
BIM.coefs$term %<>% relevel(ref="1-clone")

BIM.coefs <- BIM.coefs %>% 
  filter(term!="2 dpi", term!="3 dpi")

pd1 <- position_dodge(1)
pd2 <- position_dodge(0.7)
p3 <- ggplot(aes(y=mode, x=term), data=BIM.coefs)+
  geom_errorbar(aes(ymin=l.67, ymax=h.67), width=0, size=4, alpha=0.5)+
  geom_errorbar(aes(ymin=l.89, ymax=h.89), width=0, size=2, alpha=0.5)+
  geom_errorbar(aes(ymin=l.95, ymax=h.95), width=0)+
  geom_point(fill="white", pch=21, colour="black", size=2)+
  geom_hline(yintercept=0, linetype=2)+
  #coord_cartesian(ylim = c( -0.2251005, 1.6714590))+
  #scale_y_continuous(breaks=seq(-0.5, 2, 0.5))+
  #coord_flip()+
  cowplot::theme_cowplot()+
  labs(y="Susceptible clone selection rate", x="Diversity")+ 
  theme(axis.text = element_text(size=12),
        axis.title = element_text(face="bold", size=16))+
  NULL
p3


#### ---- Plotting predictions: BIMs ---- ####

d.BIM$m3.means <- predict(m3, d.BIM, type="response")
d.BIM$m4.means <- predict(m4, d.BIM, type="response")

d.BIM$m4.95.lower <- predict(m4, d.BIM, type="response",
                             interval="confidence", level=0.95)[,"lwr"]
d.BIM$m4.95.upper <- predict(m4, d.BIM, type="response",
                             interval="confidence", level=0.95)[,"upr"]
d.BIM$m4.89.lower <- predict(m4, d.BIM, type="response",
                             interval="confidence", level=0.89)[,"lwr"]
d.BIM$m4.89.upper <- predict(m4, d.BIM, type="response",
                             interval="confidence", level=0.89)[,"upr"]
d.BIM$m4.67.lower <- predict(m4, d.BIM, type="response",
                             interval="confidence", level=0.67)[,"lwr"]
d.BIM$m4.67.upper <- predict(m4, d.BIM, type="response",
                             interval="confidence", level=0.67)[,"upr"]

BIM.preds <- d.BIM %>% 
  select(Treatment, Timepoint, m4.means, 
         m4.95.lower, m4.95.upper,
         m4.89.lower, m4.89.upper,
         m4.67.lower, m4.67.upper) %>% 
  unique %>% 
  filter(Timepoint!="0")

pd <- position_dodge(0.5)

p4 <- ggplot(aes(x=Treatment, y=m4.means), data=BIM.preds)+
  labs(x="Diversity", y="Susceptible clone selection rate")+
  geom_hline(yintercept = 0, linetype=2, colour="black")+
  geom_errorbar(aes(ymin=m4.67.lower, ymax=m4.67.upper),
                width=0, size=4, alpha=0.5, position=pd)+
  geom_errorbar(aes(ymin=m4.89.lower, ymax=m4.89.upper),
                width=0, size=2, alpha=0.5, position=pd)+
  geom_errorbar(aes(ymin=m4.95.lower, ymax=m4.95.upper), 
                width=0, position=pd)+
  geom_point(aes(position=Treatment), pch=21, fill="white",
             colour="black", position=pd, size=2)+
  facet_wrap(~Timepoint)+
  cowplot::theme_cowplot()+
  scale_y_continuous(breaks=seq(-4,2,1))+
  coord_cartesian(ylim=c(-4.226857, 2.425297))+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(face="bold", size=16),
        #panel.grid.major = element_line(colour="lightgrey", size=.25),
        strip.text = element_text(face="bold"))+
  NULL

p4

#### ---- Save plots ---- ####
# p1 == CRISPR coefficient plot
# p2 == CRISPR predictions plot
# p3 == BIM coefficient plot
# p4 == BIM predictions plot
library(cowplot)

Fig3 <- plot_grid(p1+labs(y="Selection rate", x="")+theme(axis.text.x=element_blank()), 
                     p3+labs(y="Selection rate", x="CRISPR allele diversity"),
                     labels=c("A", "B"), ncol=1, label_size = 18)
last_plot()

ggsave("Figure_3.png", Fig3, path="~/Documents/OneDrive - University of Exeter/Papers/Common, Walker-Sunderhauf and Westra 2019/",
       device="png", dpi=600, width=18, height = 23, units=c("cm"))

FigS3 <- CRISPR_fitness_plots

ggsave("Figure_S3.png", FigS3, path="~/Documents/OneDrive - University of Exeter/Papers/Common, Walker-Sunderhauf and Westra 2019/",
       device="png",dpi=600, width=30, height = 30, units=c("cm"))

FigS4 <- BIM_fitness_plots

FigS4
ggsave("Figure_S4.png", FigS4, path="~/Documents/OneDrive - University of Exeter/Papers/Common, Walker-Sunderhauf and Westra 2019/", device="png",
       dpi=600, width=30, height = 30, units=c("cm"))

