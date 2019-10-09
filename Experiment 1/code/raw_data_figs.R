#### ---- Relative diversity Exp 1:  ---- ####
#### Created: 04/10/19 by Jack Common

rm(list=ls())

#### ---- Dependencies ---- ####
library(tidyverse)
library(magrittr)
library(scales)
library(ggstatsplot)
library(lme4)
library(cowplot)
library(lmerTest)
library(ggdark)

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
    "3" = "3-host x 3-phage",
    "6" = "6-host x 6-phage",
    "12" = "12-host x 12-phage",
    "24" = "24-host x 24-phage"
  )
  
treatment_labeller <- function(variable, value) {
    return(treatment_names[value])
}

treatment_names <- c("3-host x 3-phage", "6-host x 6-phage", "12-host x 12-phage", "24-host x 24-phage")

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
data_original <- read.csv("./Experiment 1/original_data/reldiv_exp1_master.csv", header=T, skip=1)  %>% 
  select(Treatment, Timepoint, Replicate, BIM_mix, Escape_phage, Tracked_BIM,
         PFU, SM_CFU, CRISPR_CFU, BIM_CFU, total_CFU,
         rCRISPR, rBIM)

data_original$Treatment %<>% as.factor()
data_original$Timepoint %<>% as.factor()
data_original$Replicate %<>% as.factor
data_original$Tracked_BIM %<>% as.factor()

# Host titre dataset
host <- data_original %>% 
  select(-PFU, -rCRISPR, -rBIM) %>% 
  rename(SM=SM_CFU, CRISPR=CRISPR_CFU, BIM=BIM_CFU, Total=total_CFU) %>% 
  gather("SM", "CRISPR", "BIM", "Total", key="Strain", value="Titre", factor_key = T)

# Phage titre dataset
phage <- data_original %>% 
  select(-SM_CFU, -CRISPR_CFU, -BIM_CFU, -rCRISPR, -rBIM) %>% 
  gather("PFU", key="Strain", value="Titre", factor_key = T)


# Relevel the datasets so things make sense
data_original$Treatment %<>% relevel(ref="24")
data_original$Treatment %<>% relevel(ref="12")
data_original$Treatment %<>% relevel(ref="6")
data_original$Treatment %<>% relevel(ref="3")

host$Treatment %<>% relevel(ref="24")
host$Treatment %<>% relevel(ref="12")
host$Treatment %<>% relevel(ref="6")
host$Treatment %<>% relevel(ref="3")

phage$Treatment %<>% relevel(ref="24")
phage$Treatment %<>% relevel(ref="12")
phage$Treatment %<>% relevel(ref="6")
phage$Treatment %<>% relevel(ref="3")


write.csv(phage, "./Experiment 1/original_data/phage_data.csv", row.names = F)

#### ---- Visualise data ---- ####
#### -- Phage titre -- #####
phage_plot <- ggplot(aes(y=log10(Titre+1), x=Timepoint, group=Replicate), data=phage)+
  
  #geom_point(stat='identity')+
  geom_path(stat='identity')+
  geom_hline(yintercept = 100, linetype=2)+
  facet_wrap(~Treatment, labeller = treatment_labeller)+
  
  labs(x="Days post-infection", y=expression(bold("Pfu ml"*{}^{-1}*"")))+
  #ggtitle("Density of phage")+
  
  cowplot::theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16),
        axis.title = element_text(face="bold", size=16),
        strip.text = element_text(face='bold', size=14, lineheight = 3))+
  
  coord_cartesian(ylim=c(0, 12))+
  scale_y_continuous(breaks=c(seq(0,12,2)),
                     labels=pfu_labels)+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))+
  NULL

last_plot()

#

#### -- CRISPR and SM titre -- ####
CRISPR_SM_plot <- ggplot(aes(y=Titre+1, x=Timepoint, group=Replicate), data=filter(host, Strain==c("CRISPR", "SM")))+
  
  #geom_point(stat='identity')+
  geom_line(stat='identity', aes(colour=Strain))+
  facet_wrap(~Treatment, labeller = treatment_labeller)+
  
  labs(x="Days post-infection", y=expression(bold("Cfu ml"*{}^{-1}*"")))+
  ggtitle("Density of CRISPR and SM clones ")+
  
  cowplot::theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16),
        axis.title = element_text(face="bold", size=16),
        strip.text = element_text(face='bold', size=14))+
  
  # scale_y_continuous(breaks=c(seq(0,12,1)))+
  # coord_cartesian(ylim=c(0, 12))+
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))+
  scale_colour_manual(values=c("black", "blue"))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))+
  NULL

last_plot()

#

#### -- CRISPR and BIM titre -- ####
CRISPR_BIM_plot <- ggplot(aes(y=Titre+1, x=Timepoint, group=Replicate), data=filter(host, Strain==c("CRISPR", "BIM")))+
  
  #geom_point(stat='identity')+
  geom_line(stat='identity', aes(colour=Strain))+
  facet_wrap(~Treatment, labeller = treatment_labeller)+
  
  labs(x="Days post-infection", y=expression(bold("Cfu ml"*{}^{-1}*"")))+
  ggtitle("Density of CRISPR clones and labelled BIM")+
  
  cowplot::theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16),
        axis.title = element_text(face="bold", size=16),
        strip.text = element_text(face='bold', size=14))+
  
  # scale_y_continuous(breaks=c(seq(0,12,1)))+
  # coord_cartesian(ylim=c(0, 12))+
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)))+
  scale_colour_manual(values=c("black", "blue"))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))+
  NULL

last_plot()
#
#### -- Show & save plots -- ####
quartz()
phage_plot
CRISPR_SM_plot
CRISPR_BIM_plot

ggsave("phage_plot.png", phage_plot, path="Experiment 1/docs/", device="png", dpi=600,
       width=28, height = 20, units = c("cm"))

ggsave("cfu_CRISPR_SM.png", CRISPR_SM_plot, path="./Experiment 1/docs/", device="png", dpi=600,
       width=28, height = 20, units = c("cm"))

ggsave("cfu_CRISPR_BIM.png", CRISPR_BIM_plot, path="./Experiment 1/docs/", device="png", dpi=600,
       width=28, height = 20, units = c("cm"))






