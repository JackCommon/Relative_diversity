library(data.table)
library(ggplot2)
library(plyr)
library(ape)
library(coda)
library(lattice)
library(VCVglmm)
library(MCMCglmm)
library(lme4)
library(aod)
library(magrittr)
library(tidyverse)

#### ---- Data ---- ####

data_original <- read.csv("./Experiment 1/original_data/reldiv_exp1_master.csv", header=T, skip=1)  %>% 
  select(Treatment, Timepoint, Replicate, BIM_mix, Escape_phage, Tracked_BIM,
         PFU, SM_CFU, CRISPR_CFU, BIM_CFU, white_CFU, total_CFU,
         w_SM, w_CRISPR, w_BIM, w_white)

data_original$Timepoint %<>% as.factor()
data_original$Replicate %<>% as.factor
data_original$Tracked_BIM %<>% as.factor()

# Phage titre dataset
phage <- data_original %>% 
  select(-SM_CFU, -CRISPR_CFU, -BIM_CFU, -white_CFU, -w_SM, -w_CRISPR, -w_BIM, -w_white) %>% 
  gather("PFU", key="Strain", value="Titre", factor_key = T)

phage$Treatment %<>% relevel(ref="24")
phage$Treatment %<>% relevel(ref="12")
phage$Treatment %<>% relevel(ref="6")
phage$Treatment %<>% relevel(ref="3")

phage <- phage %>% na.exclude()

#### ---- Models ---- ####

m1 <- MCMCglmm(log(PFU+1)~Treatment, random=~Replicate, data=data_original,
               nitt = 13000, burnin = 3000, thin = 10, pr = TRUE)

m2 <- MCMCglmm(log(PFU+1)~Timepoint, random=~Replicate, data=data_original,
               nitt = 13000, burnin = 3000, thin = 10, pr = TRUE)

m3 <- MCMCglmm(log(PFU+1)~Treatment + Timepoint, random=~Replicate, data=data_original,
               nitt = 13000, burnin = 3000, thin = 10, pr = TRUE)

m4 <- MCMCglmm(log(PFU+1)~Treatment * Timepoint, random=~Replicate, data=data_original,
               nitt = 13000, burnin = 3000, thin = 10, pr = TRUE)
m1$DIC
m2$DIC
m3$DIC
m4$DIC

par(mfrow=c(2,2))
MCMCfixplot(m1)
MCMCfixplot(m2)
MCMCfixplot(m3)
MCMCfixplot(m4)

summary(m3)$solutions %>% nrow()

Wald.test.auto(m1)
Wald.test.auto(m2)
Wald.test.auto(m3)
Wald.test.auto(m4)

plot(m3)
Solapply(m1); Solapply(m2); Solapply(m3)
MCMCRepnorm2(m1); MCMCRepnorm2(m2); MCMCRepnorm2(m3); MCMCRepnorm2(m4)

posterior.mode(m1$VCV); posterior.mode(m2$VCV); posterior.mode(m3$VCV)
HPDinterval(m1$VCV); HPDinterval(m2$VCV)
plot(m1$VCV)

summary(m3)

plot(m3$VCV)

#### ---- Figures ----####

# Set up target dataframe to store modes and HPDs

coefs <- data.frame(term = factor(10), 
                    mode = numeric(10), 
                    l.95 = numeric(10), h.95 = numeric(10), 
                    l.89 = numeric(10), h.89 = numeric(10), 
                    l.67 = numeric(10), h.67 = numeric(10))

coefs$term <- c("Intercept", "6", "12", "24","1 dpi", "2 dpi", "3 dpi") %>% 
  as.factor

# Get posterior coefficients
coefs$mode <- posterior.mode(m3$Sol)[1:10]
coefs$l.95 <- HPDinterval(m3$Sol, prob=0.95) %>% 
  as.data.frame() %>%    
  slice(1:10) %>% 
  select(lower) %>% 
  .[1:10,1]
coefs$h.95 <- HPDinterval(m3$Sol, prob=0.95) %>% 
  as.data.frame() %>%    
  slice(1:10) %>% 
  select(upper) %>% 
  .[1:10,1]

coefs$l.89 <- HPDinterval(m3$Sol, prob=0.89) %>% 
  as.data.frame() %>%    
  slice(1:10) %>% 
  select(lower ) %>% 
  .[1:10,1]
coefs$h.89 <- HPDinterval(m3$Sol, prob=0.89) %>% 
  as.data.frame() %>%    
  slice(1:10) %>% 
  select(upper) %>% 
  .[1:10,1]

coefs$l.67 <- HPDinterval(m3$Sol, prob=0.67) %>% 
  as.data.frame() %>%    
  slice(1:10) %>% 
  select(lower) %>% 
  .[1:10,1]
coefs$h.67 <- HPDinterval(m3$Sol, prob=0.67) %>% 
  as.data.frame() %>%    
  slice(1:10) %>% 
  select(upper) %>% 
  .[1:10,1]

coefs$term %<>% relevel(ref="3 dpi")
coefs$term %<>% relevel(ref="2 dpi")
coefs$term %<>% relevel(ref="1 dpi")
coefs$term %<>% relevel(ref="24")
coefs$term %<>% relevel(ref="12")
coefs$term %<>% relevel(ref="6")
coefs$term %<>% relevel(ref="Intercept")

pd1 <- position_dodge(1)
pd2 <- position_dodge(0.7)
p1 <- ggplot(aes(y=mode, x=term), data=coefs)+
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
p1

ggsave("phage_model_coefs.png", p1, path='./Experiment 1/docs/', device="png", dpi=600,
       width=20, height=12, units=c("cm"))

#ggsave("Figure_S2.png", p1, path='~/Documents/OneDrive - University of Exeter/Papers/Common, Walker-Sunderhauf and Westra 2019/', 
#       device="png", dpi=600,width=20, height=12, units=c("cm"))

