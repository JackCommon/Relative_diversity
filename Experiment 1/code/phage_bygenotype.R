#### phage_bygenotype.R by Jack Common
#### Code to generate figures that display the raw phage titres of each individual phage genotype
#### in each diversity treatment. 

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


#### ---- Data ---- ####
phage <- read.csv("./Experiment 1/original_data/phage_bygenotype_data.csv", header=T, skip=1)

phage$Treatment %<>% as.character()
phage$Replicate %<>% as.character()
phage$Replicate_new <- paste(phage$Treatment, ".", phage$Replicate, ".", phage$Genotype)

phage$Timepoint %<>% as.factor()
phage$Genotype %<>% as.factor()
phage$Treatment %<>% as.factor()
phage$Replicate_new %<>% as.factor()

phage$Treatment %<>% relevel(ref="24")
phage$Treatment %<>% relevel(ref="12")
phage$Treatment %<>% relevel(ref="6")
phage$Treatment %<>% relevel(ref="3")

phage %<>% na.exclude()

#### ---- Figure: facetted by treatment ---- ####
# First plot the titres of each genotype
bygenotype_plot <- ggplot(aes(y=log10(PFU+1), x=Timepoint), data=phage)+
  
  #geom_point(stat='identity')+
  geom_hline(yintercept = 1, linetype=2)+
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
  scale_colour_discrete(name=c("Phage genotype"))+
  
  theme(axis.text = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16, face="bold"))+
  NULL

last_plot() 

### ---- now let's bring in the overall phage titre data... ---- ####
phage_all <- read.csv("./Experiment 1/original_data/phage_data.csv")

phage_all$Timepoint %<>% as.factor()
phage_all$Treatment %<>% as.factor()
phage_all$Replicate %<>% as.factor()

phage_all$Treatment %<>% relevel(ref="24")
phage_all$Treatment %<>% relevel(ref="12")
phage_all$Treatment %<>% relevel(ref="6")
phage_all$Treatment %<>% relevel(ref="3")

phage_all <- phage_all %>% group_by(Treatment, Timepoint) %>% 
  mutate(mid=median(Titre), min=min(Titre), max=max(Titre))

phage_all <- phage_all %>% rename(PFU=Titre)

bygenotype_plot_2 <- bygenotype_plot + 
  geom_ribbon(data=filter(phage_all, Timepoint!="0"), fill="grey", alpha=0.5,
              aes(ymin=log10(min), ymax=log10(max), group=Treatment))+
  geom_path(stat='identity', aes(colour=Genotype, group=Replicate_new))+
  geom_path(stat='identity', aes(y=log10(mid), x=Timepoint, group=Treatment),
            data=filter(phage_all, Timepoint!="0"), size=1)+
  NULL

last_plot()

ggsave("phage_bygenotype.png", bygenotype_plot_2, path="./Experiment 1/figs/",
       device="png", dpi=600, width=28, height=20, units=c("cm"))

#### ---- Figure: look at each replicate in each treatment ---- ####

# Make subset dataframes by treatment
three <- phage %>% 
  filter(Treatment=="3")
six <- phage %>% 
  filter(Treatment=="6")
twelve <- phage %>% 
  filter(Treatment=="12")
twentyfour <- phage %>% 
  filter(Treatment=="24")

# Set up blank plots for each treatment
three_bygenotype_plot <- ggplot(aes(y=log10(PFU+1), x=Timepoint), data=three)+
  
  #geom_point(stat='identity')+
  geom_hline(yintercept = 1, linetype=2)+
  facet_wrap(~Replicate, nrow=2, ncol=4)+
  
  labs(x="Days post-infection", y=expression(bold("Pfu ml"*{}^{-1}*"")))+
  ggtitle("A. 3-clone x 3-phage")+
  
  cowplot::theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16),
        axis.title = element_text(face="bold", size=16),
        strip.text = element_text(face='bold', size=14, lineheight = 3))+
  
  coord_cartesian(ylim=c(0, 12))+
  scale_y_continuous(breaks=c(seq(0,12,2)),
                     labels=pfu_labels)+
  scale_colour_discrete(name=c("Phage genotype"))+
  
  theme(axis.text = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16, face="bold"))+
  NULL

six_bygenotype_plot <- ggplot(aes(y=log10(PFU+1), x=Timepoint), data=six)+
  
  #geom_point(stat='identity')+
  geom_hline(yintercept = 1, linetype=2)+
  facet_wrap(~Replicate, nrow=2, ncol=4)+
  
  labs(x="Days post-infection", y=expression(bold("Pfu ml"*{}^{-1}*"")))+
  ggtitle("B. 6-clone x 6-phage")+
  
  cowplot::theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16),
        axis.title = element_text(face="bold", size=16),
        strip.text = element_text(face='bold', size=14, lineheight = 3))+
  
  coord_cartesian(ylim=c(0, 12))+
  scale_y_continuous(breaks=c(seq(0,12,2)),
                     labels=pfu_labels)+
  scale_colour_discrete(name=c("Phage genotype"))+
  
  theme(axis.text = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16, face="bold"))+
  NULL

twelve_bygenotype_plot <- ggplot(aes(y=log10(PFU+1), x=Timepoint), data=twelve)+
  
  #geom_point(stat='identity')+
  geom_hline(yintercept = 1, linetype=2)+
  facet_wrap(~Replicate, nrow=2, ncol=4)+
  
  labs(x="Days post-infection", y=expression(bold("Pfu ml"*{}^{-1}*"")))+
  ggtitle("C. 12-clone x 12-phage")+
  
  cowplot::theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16),
        axis.title = element_text(face="bold", size=16),
        strip.text = element_text(face='bold', size=14, lineheight = 3))+
  
  coord_cartesian(ylim=c(0, 12))+
  scale_y_continuous(breaks=c(seq(0,12,2)),
                     labels=pfu_labels)+
  scale_colour_discrete(name=c("Phage genotype"))+
  
  theme(axis.text = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16, face="bold"))+
  NULL

twentyfour_bygenotype_plot <- ggplot(aes(y=log10(PFU+1), x=Timepoint), data=twentyfour)+
  
  #geom_point(stat='identity')+
  geom_hline(yintercept = 1, linetype=2)+
  facet_wrap(~Replicate, nrow=2, ncol=4)+
  
  labs(x="Days post-infection", y=expression(bold("Pfu ml"*{}^{-1}*"")))+
  ggtitle("D. 24-clone x 24-phage")+
  
  cowplot::theme_cowplot()+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16),
        axis.title = element_text(face="bold", size=16),
        strip.text = element_text(face='bold', size=14, lineheight = 3))+
  
  coord_cartesian(ylim=c(0, 12))+
  scale_y_continuous(breaks=c(seq(0,12,2)),
                     labels=pfu_labels)+
  scale_colour_discrete(name=c("Phage genotype"))+
  
  theme(axis.text = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16, face="bold"))+
  NULL

#

three_bygenotype_plot <- three_bygenotype_plot + 
  geom_path(stat='identity', aes(colour=Genotype, group=Replicate_new))+
  theme(legend.position = "none")+
  NULL
last_plot()

six_bygenotype_plot <- six_bygenotype_plot + 
  geom_path(stat='identity', aes(colour=Genotype, group=Replicate_new))+
  theme(legend.position = "none")+
  NULL
last_plot()

twelve_bygenotype_plot <- twelve_bygenotype_plot + 
  geom_path(stat='identity', aes(colour=Genotype, group=Replicate_new))+
  theme(legend.position = "none")+
  NULL
last_plot()

twentyfour_bygenotype_plot <- twentyfour_bygenotype_plot + 
  geom_path(stat='identity', aes(colour=Genotype, group=Replicate_new))+
  theme(legend.position = "none")+
  NULL
last_plot()

byreplicate_plot <- plot_grid(three_bygenotype_plot+xlab(""),
                              six_bygenotype_plot+labs(y="", x=""),
                              twelve_bygenotype_plot, twentyfour_bygenotype_plot+ylab(""))
quartz()
last_plot()

ggsave("phage_bygenotype_byreplicate.png", byreplicate_plot, path="./Experiment 1/figs/",
       device="png", dpi=600, width=50, height=30, units=c("cm"))
