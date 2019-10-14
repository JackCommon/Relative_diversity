#### phage_analysis.R by Jack Common
#### Code to generate models that compare phage titre as a function of days post-infection (dpi) and
#### phage x host diversity treatment

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
phage <- read.csv("Experiment 1/original_data/phage_data.csv", header=T) %>% 
  select(-Escape_phage)
phage$Treatment %<>% as.factor
phage$Timepoint %<>% as.factor

phage$Treatment %<>% relevel(ref="24")
phage$Treatment %<>% relevel(ref="12")
phage$Treatment %<>% relevel(ref="6")
phage$Treatment %<>% relevel(ref="3")

phage$Timepoint %<>% relevel(ref="3")
phage$Timepoint %<>% relevel(ref="2")
phage$Timepoint %<>% relevel(ref="1")
phage$Timepoint %<>% relevel(ref="0")

#### --- lme4 models ---- ####
m1 <- lmer(log(Titre+1)~Treatment+(1|Replicate), data=phage)
m2 <- lmer(log(Titre+1)~Timepoint+(1|Replicate), data=phage)
m3 <- lmer(log(Titre+1)~Treatment+Timepoint+(1|Replicate), data=phage)
m4 <- lmer(log(Titre+1)~Treatment*Timepoint+(1|Replicate), data=phage)

AIC(m1, m2, m3, m4) %>% compare_AICs()
plot(m1)
plot(m2)
plot(m3)
plot(m4)

anova(m1, m2, m3, m4, type="marginal")
anova(m3, type="marginal")
summary(m3)

#### ---- Coefficient plot ---- ####
m3 <- lmer(log(Titre+1)~Treatment + Timepoint+(1|Replicate), data=phage)

coefs <- data.frame(term = factor(7), 
                    beta = numeric(7), 
                    l.95 = numeric(7), h.95 = numeric(7), 
                    l.89 = numeric(7), h.89 = numeric(7), 
                    l.67 = numeric(7), h.67 = numeric(7))

coefs$term <- c("Intercept", "6", "12", "24", "1 dpi", "2 dpi", "3 dpi") %>% 
  as.factor

# Get coefficients and confidence intervals
coefs$beta <- fixef(m3)
coefs$l.95 <- confint(m3, parm="beta_", level=0.95) %>% 
  as.data.frame() %>%    
  slice(1:7) %>% 
  select("2.5 %") %>% 
  .[1:7,1]
coefs$h.95 <- confint(m3, parm="beta_", level=0.95) %>% 
  as.data.frame() %>%    
  slice(1:7) %>% 
  select("97.5 %") %>% 
  .[1:7,1]

coefs$l.89 <- confint.merMod(m3, parm="beta_", level=0.89) %>% 
  as.data.frame() %>%    
  slice(1:7) %>% 
  select("5.5 %") %>% 
  .[1:7,1]
coefs$h.89 <- confint(m3, parm="beta_", level=0.89) %>% 
  as.data.frame() %>%    
  slice(1:7) %>% 
  select("94.5 %") %>% 
  .[1:7,1]

coefs$l.67 <- confint.merMod(m3, parm="beta_", level=0.67) %>% 
  as.data.frame() %>%    
  slice(1:7) %>% 
  select("16.5 %") %>% 
  .[1:7,1]
coefs$h.67 <- confint(m3, parm="beta_", level=0.67) %>% 
  as.data.frame() %>%    
  slice(1:7) %>% 
  select("83.5 %") %>% 
  .[1:7,1]

write.csv(coefs, "./Experiment 1/summary_data/phage_model_coefs.csv", row.names = F)

coefs$term %<>% relevel(ref="3 dpi")
coefs$term %<>% relevel(ref="2 dpi")
coefs$term %<>% relevel(ref="1 dpi")
coefs$term %<>% relevel(ref="24")
coefs$term %<>% relevel(ref="12")
coefs$term %<>% relevel(ref="6")
coefs$term %<>% relevel(ref="Intercept")

phage_coefs <- ggplot(aes(y=beta, x=term), data=coefs)+
  geom_errorbar(aes(ymin=l.67, ymax=h.67), width=0, size=4, alpha=0.5)+
  geom_errorbar(aes(ymin=l.89, ymax=h.89), width=0, size=2, alpha=0.5)+
  geom_errorbar(aes(ymin=l.95, ymax=h.95), width=0)+
  geom_point(fill="white", pch=21, colour="black", size=2)+
  geom_hline(yintercept=0, linetype=2)+
  coord_flip()+
  cowplot::theme_cowplot()+
  labs(y=expression(bold("ln(pfu ml"*{}^{-1}*")")), x="Fixed effect level")+ 
  scale_y_continuous(breaks=seq(-14, 20, 2))+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(face="bold", size=16))+
  NULL
last_plot()

ggsave("phage_coefs.png", p1, path="./Experiment 1/figs/", 
       device="png", dpi=600,width=15, height=12, units=c("cm"))
