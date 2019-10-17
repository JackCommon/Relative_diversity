#### comparison_code.R by Jack Common
#### Code to generate figures and models that compare phage titre and selection rates
#### between experiments where a gradient of polyclonal CRISPR host populations were  
#### treated with either monoclonal or polyclonal phage. Data for monoclonal phage
#### experiments is derived from Common et al 2019 Ecol Lett

rm(list=ls())
options(show.signif.stars = F)

#### ---- Dependencies ---- ####
library(scales)
library(cowplot)
library(ggplot2)
library(lme4)
library(magrittr)
library(tidyverse)
library(lmerTest)

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
  "3" = "3",
  "6" = "6",
  "12" = "12",
  "24" = "24"
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

treatment_list <- c("3", "6", "12", "24")

pfu_labels <- c(
  expression("0"), 
  expression('10'^2*''), 
  expression('10'^4*''),
  expression('10'^6*''),
  expression('10'^8*''),
  expression('10'^10*''),
  expression("10"^12*"")
)

pd <- position_dodge(0.5)

##### ---- Phage data ---- #####
reldiv_phage <- read.csv("Experiment 1/original_data/phage_data.csv", header=T) %>% 
  select(-Escape_phage)
reldiv_phage$Experiment <- rep("reldiv", nrow(reldiv_phage))
reldiv_phage$Treatment %<>% as.factor

ecolet_phage <- read.csv("Experiment 1/original_data/ecolet_phage_data.csv", header=T) %>% 
  select(-Escape_phage) %>% 
  filter(Treatment %in% c("3", "6", "12", "24"))
ecolet_phage$Experiment <- rep("ecolet", nrow(ecolet_phage))
ecolet_phage$Treatment %<>% as.factor()

combined_phage <- bind_rows(reldiv_phage, ecolet_phage)
combined_phage$Treatment %<>% as.factor
combined_phage$Timepoint %<>% as.factor
combined_phage$Replicate %<>% as.factor
combined_phage$Experiment %<>% as.factor

# Relevel the datasets so things make sense
combined_phage$Treatment %<>% relevel(ref="24")
combined_phage$Treatment %<>% relevel(ref="12")
combined_phage$Treatment %<>% relevel(ref="6")
combined_phage$Treatment %<>% relevel(ref="3")

#### -- Phage titre plot-- #####
phage_comparison <- ggplot(aes(y=log10(Titre+1), x=Timepoint), data=combined_phage)+
  
  geom_boxplot(width=0.5, outlier.shape = NA, position=pd, aes(colour=Experiment))+
  geom_point(width=0.2, alpha=0.5, pch=21, position=pd, aes(colour=Experiment))+
  # stat_summary(fun.y="mean", geom="point", aes(colour=Experiment), position=pd)+
  # stat_summary(fun.data = "mean_cl_boot", geom="errorbar", width=0, aes(colour=Experiment), position=pd)+
  #geom_path(stat='identity', alpha=0.5, position=pd, aes(colour=Experiment, group=Replicate))+
  #geom_hline(yintercept = 100, linetype=2)+
  facet_wrap(~Treatment)+
  
  labs(x="Days post-infection", y=expression(bold("Pfu ml"*{}^{-1}*"")))+
  #ggtitle("Density of phage")+
  
  cowplot::theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16),
        axis.title = element_text(face="bold", size=16),
        strip.text = element_text(face='bold', size=14, lineheight = 3))+
  
  coord_cartesian(ylim=c(0, 12))+
  scale_y_continuous(breaks=c(seq(0,12,2)),
                     labels=pfu_labels)+
  scale_colour_discrete(labels = c("CRISPR allele diversity", "CRISPR allele & phage\ndiversity"))+
  
  theme(axis.text = element_text(size=12),
        legend.title = element_text(face="bold", size=16),
        legend.text = element_text(size=12),
        legend.key.height = unit(1, "cm"),
        panel.grid.major = element_line(colour="lightgrey", size=.25))+
  NULL

last_plot()

ggsave("phage_comparison.png", phage_comparison, path="Experiment 1/figs/", device="png", dpi=600,
       width=28, height = 20, units = c("cm"))

#

#### ---- Phage titre lme4 ---- ####
m1 <- lmer(log(Titre+1)~Experiment+Treatment+Timepoint+(1|Replicate), 
           data=filter(combined_phage, Timepoint!="0"))
summary(m1)
anova(m1)

confint(m1, parm="beta_", method="boot")

#### ---- Selection rates ---- ####
reldiv_rates <- read.csv("./Experiment 1/original_data/reldiv_exp1_master.csv", header=T, skip = 1)  %>% 
  select(Treatment, Timepoint, Replicate, BIM_mix, Escape_phage, Tracked_BIM,
         rCRISPR, rBIM) %>% 
  gather("rCRISPR", "rBIM", key="Strain", value="r", factor_key = T) %>% 
  filter(Timepoint %in% c("1", "2", "3"), r!='#NUM!')

reldiv_rates$Experiment <- rep("reldiv", nrow(reldiv_rates))

reldiv_rates$Timepoint %<>% as.factor()
reldiv_rates$Replicate %<>% as.factor
reldiv_rates$Tracked_BIM %<>% as.factor()
reldiv_rates$Escape_phage %<>% as.factor
reldiv_rates$Treatment %<>% as.factor
reldiv_rates$r %<>% as.numeric()

ecolet_rates <- read.csv("./Experiment 1/original_data/ecolet_rates.csv", header=T)  %>% 
  select(Treatment, Timepoint, Replicate, BIM_mix, Escape_phage, Tracked_BIM,
         rCRISPR, rBIM) %>% 
  gather("rCRISPR", "rBIM", key="Strain", value="r", factor_key = T) %>% 
  filter(Timepoint %in% c("1", "2", "3"), Treatment %in% c("3", "6", "12", "24"))

ecolet_rates$Timepoint %<>% as.factor()
ecolet_rates$Replicate %<>% as.factor
ecolet_rates$Tracked_BIM %<>% as.factor()
ecolet_rates$Escape_phage %<>% as.factor
ecolet_rates$Treatment %<>% as.factor

ecolet_rates$Experiment <- rep("ecolet", nrow(ecolet_rates))

combined_rates <- bind_rows(reldiv_rates, ecolet_rates)
combined_rates$Treatment %<>% as.factor

# Relevel the datasets so things make sense
combined_rates$Treatment %<>% relevel(ref="24")
combined_rates$Treatment %<>% relevel(ref="12")
combined_rates$Treatment %<>% relevel(ref="6")
combined_rates$Treatment %<>% relevel(ref="3")

#### ---- CRISPR boxplot ---- ####

pd <- position_dodge(0.2)

CRISPR_rate_comp <- ggplot(aes(x=Treatment, y=r), 
                            data=filter(combined_rates, Strain=="rCRISPR", Timepoint!="0"))+
  geom_boxplot(width=0.5, outlier.shape = NA, aes(colour=Experiment))+
  geom_jitter(width=0.2, alpha=0.5, pch=21, aes(colour=Experiment))+
  geom_hline(yintercept=0, linetype=2, size=1)+
  # stat_summary(fun.y="mean", geom="point", aes(colour=Experiment, position=Experiment), position=pd)+
  # stat_summary(fun.data = "mean_cl_boot", geom="errorbar", width=0, aes(colour=Experiment), position=pd)+
   facet_wrap(~Timepoint, labeller = timepoint_labeller, ncol=2)+
  labs(y="Selection rate (day -1) of CRISPR", x="Number of CRISPR and phage clones")+
  scale_x_discrete(labels = treatment_list)+
  #scale_y_continuous(breaks=c(seq(-1.5, 2, 0.5)))+
  scale_colour_discrete(labels = c("CRISPR allele diversity", "CRISPR allele & phage\ndiversity"))+
  theme_cowplot()+
  theme(axis.title = element_text(face="bold", size=16),
        legend.title = element_text(face="bold", size=16),
        panel.grid.major = element_line(colour="lightgrey", size=.25),
        strip.text = element_text(face="bold"),
        axis.text.x = element_text(size=10),
        legend.key.height = unit(1, "cm"))+
  NULL

last_plot()

ggsave("rate_CRISPR_comparison.png", CRISPR_rate_comp, path="Experiment 1/figs/", device="png", dpi=600,
       width=28, height = 20, units = c("cm"))

#### ---- CRISPR lme4 ---- #####
d.CRISPR <- combined_rates %>% 
  na.exclude %>% 
  filter(Strain=="rCRISPR", Timepoint!="0")

m1 <- lmer(r~Experiment+(1|Replicate), 
           data=d.CRISPR)
m2 <- lmer(r~Experiment+Treatment+(1|Replicate), 
           data=d.CRISPR)
m3 <- lmer(r~Experiment+Timepoint+(1|Replicate), 
           data=d.CRISPR)
m4 <- lmer(r~Experiment+Treatment+Timepoint+(1|Replicate), 
                 data=d.CRISPR)

AIC(m1, m2, m3, m4) %>% compare_AICs()

summary(m2)
anova(m4)

confint(m2, parm="beta_", method="boot")

#### ---- BIM boxplot ---- ####
BIM_rate_comp <- ggplot(aes(x=Treatment, y=r), 
                            data=filter(combined_rates, Strain=="rBIM", Timepoint!="0"))+
  geom_boxplot(width=0.5, outlier.shape = NA, aes(colour=Experiment))+
  geom_jitter(width=0.2, alpha=0.5, pch=21, aes(colour=Experiment))+
  geom_hline(yintercept=0, linetype=2, size=1)+
  # stat_summary(fun.y="mean", geom="point", aes(colour=Experiment), position=pd)+
  # stat_summary(fun.data = "mean_cl_boot", geom="errorbar", width=0, aes(colour=Experiment), position=pd)+
  facet_wrap(~Timepoint, labeller = timepoint_labeller, ncol=2)+
  labs(y="Selection rate (day -1) of\nsusceptible clone", x="Number of CRISPR and phage clones")+
  scale_x_discrete(labels = treatment_list)+
  #scale_y_continuous(breaks=c(seq(-1.5, 2, 0.5)))+
  scale_colour_discrete(labels = c("CRISPR allele diversity", "CRISPR allele & phage\ndiversity"))+
  theme_cowplot()+
  theme(axis.title = element_text(face="bold", size=16),
        legend.title = element_text(face="bold", size=16),
        panel.grid.major = element_line(colour="lightgrey", size=.25),
        strip.text = element_text(face="bold"),
        axis.text.x = element_text(size=10),
        legend.key.height = unit(1, "cm"))+
  NULL
last_plot()

ggsave("rate_BIM_comparison.png", BIM_rate_comp, path="Experiment 1/figs/", device="png", dpi=600,
       width=28, height = 20, units = c("cm"))

#### ---- BIM lme4 ---- ####
d.BIM <- combined_rates %>% 
  na.exclude %>% 
  filter(Strain=="rBIM", Timepoint!="0")

m1 <- lmer(r~Experiment+(1|Replicate), 
           data=d.BIM)
m2 <- lmer(r~Experiment+Treatment+(1|Replicate), 
           data=d.BIM)
m3 <- lmer(r~Experiment+Timepoint+(1|Replicate), 
           data=d.BIM)
m4 <- lmer(r~Experiment+Treatment+Timepoint+(1|Replicate), 
           data=d.BIM)

AIC(m1, m2, m3, m4) %>% compare_AICs()

summary(m3)
anova(m2, type="marginal")

confint(m3, parm="beta_", method="boot")

#### ---- Paper figures ---- ####
## Figure 1 - Phage comparison 
Fig1 <- ggplot(aes(y=log10(Titre+1), x=Timepoint), data=combined_phage)+
  
  geom_boxplot(width=0.5, outlier.shape = NA, position=pd, aes(colour=Experiment))+
  geom_point(width=0.2, alpha=0.5, pch=21, position=pd, aes(colour=Experiment))+
  # stat_summary(fun.y="mean", geom="point", aes(colour=Experiment), position=pd)+
  # stat_summary(fun.data = "mean_cl_boot", geom="errorbar", width=0, aes(colour=Experiment), position=pd)+
  #geom_path(stat='identity', alpha=0.5, position=pd, aes(colour=Experiment, group=Replicate))+
  #geom_hline(yintercept = 100, linetype=2)+
  facet_wrap(~Treatment)+
  
  labs(x="Days post-infection", y=expression(bold("Pfu ml"*{}^{-1}*"")))+
  #ggtitle("Density of phage")+
  
  cowplot::theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16),
        axis.title = element_text(face="bold", size=16),
        strip.text = element_text(face='bold', size=14, lineheight = 3))+
  
  coord_cartesian(ylim=c(0, 12))+
  scale_y_continuous(breaks=c(seq(0,12,2)),
                     labels=pfu_labels)+
  scale_colour_discrete(labels = c("CRISPR allele diversity", "CRISPR allele & phage\ndiversity"))+
  
  theme(axis.text = element_text(size=12),
        legend.title = element_text(face="bold", size=16),
        legend.text = element_text(size=12),
        legend.key.height = unit(1, "cm"),
        panel.grid.major = element_line(colour="lightgrey", size=.25))+
  NULL

last_plot()

ggsave("Figure 1.png", Fig1, path="Experiment 1/figs/", device="png", dpi=600,
       width=28, height = 20, units = c("cm"))

## Figure 2 - CRIPSR rate comparison

pd <- position_dodge(0.2)

Fig2 <- ggplot(aes(x=Treatment, y=r), 
                           data=filter(combined_rates, Strain=="rCRISPR", Timepoint!="0"))+
  geom_boxplot(width=0.5, outlier.shape = NA, aes(colour=Experiment))+
  geom_point(width=0.2, alpha=0.5, pch=21, position=pd, aes(colour=Experiment))+
  # geom_jitter(width=0.2, alpha=0.5, pch=21, aes(colour=Experiment))+
  geom_hline(yintercept=0, linetype=2, size=1)+
  # stat_summary(fun.y="mean", geom="point", aes(colour=Experiment, position=Experiment), position=pd)+
  # stat_summary(fun.data = "mean_cl_boot", geom="errorbar", width=0, aes(colour=Experiment), position=pd)+
  # facet_wrap(~Timepoint, labeller = timepoint_labeller, ncol=2)+
  labs(y="Selection rate of CRISPR", x="Number of CRISPR and phage clones")+
  scale_x_discrete(labels = treatment_list)+
  #scale_y_continuous(breaks=c(seq(-1.5, 2, 0.5)))+
  scale_colour_discrete(labels = c("CRISPR allele diversity", "CRISPR allele & phage\ndiversity"))+
  theme_cowplot()+
  theme(axis.title = element_text(face="bold", size=16),
        legend.title = element_text(face="bold", size=16),
        panel.grid.major = element_line(colour="lightgrey", size=.25),
        strip.text = element_text(face="bold"),
        axis.text.x = element_text(size=10),
        legend.key.height = unit(1, "cm"))+
  NULL

last_plot()

ggsave("Figure 2.png", Fig2, path="Experiment 1/figs/", device="png", dpi=600,
       width=20, height = 10, units = c("cm"))

## Figure 3 - BIM rate comparison
Fig3 <- ggplot(aes(x=Treatment, y=r), 
                        data=filter(combined_rates, Strain=="rBIM", Timepoint!="0"))+
  geom_boxplot(width=0.5, outlier.shape = NA, aes(colour=Experiment))+
  geom_point(width=0.2, alpha=0.5, pch=21, position=pd, aes(colour=Experiment))+
  # geom_jitter(width=0.2, alpha=0.5, pch=21, aes(colour=Experiment))+
  geom_hline(yintercept=0, linetype=2, size=1)+
  # stat_summary(fun.y="mean", geom="point", aes(colour=Experiment), position=pd)+
  # stat_summary(fun.data = "mean_cl_boot", geom="errorbar", width=0, aes(colour=Experiment), position=pd)+
  # facet_wrap(~Timepoint, labeller = timepoint_labeller, ncol=2)+
  labs(y="Selection rate of\nsusceptible clone", x="Number of CRISPR and phage clones")+
  scale_x_discrete(labels = treatment_list)+
  #scale_y_continuous(breaks=c(seq(-1.5, 2, 0.5)))+
  scale_colour_discrete(labels = c("CRISPR allele diversity", "CRISPR allele & phage\ndiversity"))+
  theme_cowplot()+
  theme(axis.title = element_text(face="bold", size=16),
        legend.title = element_text(face="bold", size=16),
        panel.grid.major = element_line(colour="lightgrey", size=.25),
        strip.text = element_text(face="bold"),
        axis.text.x = element_text(size=10),
        legend.key.height = unit(1, "cm"))+
  NULL
last_plot()

ggsave("Figure 3.png", Fig3, path="Experiment 1/figs/", device="png", dpi=600,
       width=20, height = 10, units = c("cm"))
