rm(list=ls())

library(tidyverse)
library(ggplot2)
library(effects)
library(ggplot2)
library(glmmTMB)
library(viridis)
library(viridisLite)
library(MetBrewer)
library(dplyr)
library(vegan)
library(GCalignR) 
library(writexl)
library(ggplot2)
library(MetBrewer)
library(glmmTMB) 
library(ggpubr)
library(parameters)
library(emmeans)

#################################################################################################
### Objective 1.  The relative role of nutrients and defensive metabolites in bat preference ###
#################################################################################################

setwd("~/Desktop/piperine/Preference experiments")
batdata <- read.csv("preference.csv")
head(batdata)
hist(batdata$percentagefoodeaten)
batdata$food_gamma <- (batdata$percentagefoodeaten/100) + 0.0001
batdata

### Filter the data per treatment. Three different two choice treatments were conducted ###
nutrients <- batdata %>% filter(treatment %in% c("hn", "ln"))
summary(nutrients)
amides <- batdata %>% filter(treatment %in% c("ha", "la"))
summary(amides)
nutrientamides <- batdata %>% filter(treatment %in% c("ln+lamides", "hn+hamides"))
summary(nutrientamides)

#############
### GLMM ###
############

### NUTRIENTS ###
nutrients_glmm <- glmmTMB(food_gamma ~ treatment, data = nutrients, family = beta_family(link="logit"))
summary(nutrients_glmm)
shapiro.test(resid(nutrients_glmm)) # 0.1115 normal 
diagnose(nutrients_glmm) # model looks OK!
plot(allEffects(nutrients_glmm))
parameters(nutrients_glmm)
nutrients$prediction <- predict(nutrients_glmm, nutrients, type="response")
head(nutrients)
nutrients_emmeans <-emmeans(nutrients_glmm,~treatment, type="response")
nutrients_emmeans
nutrients_emmeans <- as.data.frame(nutrients_emmeans)

# Effect size nutrients, high treatment = 0.515, low = 0.237
x <- 0.515 - 0.237
effect_size_n <- (x*100)/0.237
effect_size_n # Percentage: bats had 117.2996% more in high nutrient treatment 
xx <- 0.515/0.237 
xx # Fold change: bats consume 2.17 more food in the high nutrient treatment  

### AMIDES ###
hist(amides$percentagefoodeaten)
amides_glmm <- glmmTMB(food_gamma ~ treatment, data = amides, family = beta_family(link="logit"))
summary(amides_glmm)
shapiro.test(resid(amides_glmm)) # 0.09552 normal
diagnose(amides_glmm) # model looks OK!
plot(allEffects(amides_glmm))
parameters(amides_glmm)
amides$prediction <- predict(amides_glmm, amides, type="response")
amides_emmeans <-emmeans(amides_glmm,~treatment, type="response")
amides_emmeans
amides_emmeans <- as.data.frame(amides_emmeans)

# Effect size amides, high treatment = 0.161, low = 0.564
x <- 0.564 - 0.161
effect_size_a <- (x*100)/0.161
effect_size_a # Percentage: bats had 250.3106% more in low amide treatment 
xx <- 0.564/0.161
xx # Fold change: bats consume 3.50 more food in the low amide treatment

### NUTRIENTS AND AMIDES ###
hist(nutrientamides$percentagefoodeaten)
nutrientamides_glmm <- glmmTMB(food_gamma ~ treatment, data = nutrientamides, family = beta_family(link="logit"))
summary(nutrientamides_glmm)
shapiro.test(resid(nutrientamides_glmm)) #  0.02216 normal
diagnose(nutrientamides_glmm) # model looks OK!
plot(allEffects(nutrientamides_glmm))
parameters(nutrientamides_glmm)
nutrientamides$prediction <- predict(nutrientamides_glmm, nutrientamides, re.form=NA,type="response")
nutrientamides_emmeans <-emmeans(nutrientamides_glmm,~treatment, type="response")
nutrientamides_emmeans
nutrientamides_emmeans <- as.data.frame(nutrientamides_emmeans)

# Effect size amides, high N+A treatment = 0.552, N+A low = 0.172
x <- 0.552 - 0.172
effect_size_na <- (x*100)/0.172
effect_size_na # Percentage: bats had 220.93% more in the high N+A treatment 
xx <- 0.552/0.172
xx # Fold change: bats consume 3.20 more food in in the high N+A treatment 

##############
### Graphs ###
##############

nutrients_graph <- ggplot(nutrients, aes(x = treatment, y = prediction, color= treatment)) + 
  theme_classic(base_size = 13) +
  geom_jitter(data = nutrients, aes(x = treatment, y = food_gamma, color = treatment), width = 0.1, size = 3) +
  scale_color_viridis(option = "D", discrete=TRUE, direction = -1) +
  geom_errorbar(data = nutrients_emmeans, aes(x = treatment, y = response, ymin = lower.CL, ymax = upper.CL), width = 0.2, color = "black") + 
  stat_summary(fun.data = mean_se, color = "black") +
  ylab (" ") +
  xlab (" ")

a <- nutrients_graph +
  theme(legend.position = "none") + 
  scale_x_discrete(labels = c("hn" = "High nutrients", "ln" = "Low nutrients"),
                   limits = c("ln","hn"))
a

amides_graph <- ggplot(amides, aes(x = treatment, y = prediction, color= treatment)) +
  theme_classic(base_size = 13) +
  geom_jitter(data = amides, aes(x = treatment, y = food_gamma, color = treatment), width = 0.1, size = 3) +
  scale_color_viridis(option = "D", discrete=TRUE, direction = -1) +
  geom_errorbar(data = amides_emmeans, aes(x = treatment, y = response, ymin = lower.CL, ymax = upper.CL), width = 0.2, color = "black") + 
  stat_summary(fun.data = mean_se, color = "black") +
  ylab ("Proportion of food eaten") +
  xlab (" ")
b <- amides_graph +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("ha" = "2% piperine", "la" = "0.1% piperine"),
                   limits = c("la","ha"))
b

nutrientamides_graph <- ggplot(nutrientamides, aes(x = treatment, y = prediction, color= treatment)) +
  theme_classic(base_size = 13) +
  geom_jitter(data = nutrientamides, aes(x = treatment, y = food_gamma, color = treatment), width = 0.1, size = 3) +
  scale_color_viridis(option = "D", discrete=TRUE, direction = -1) +
  geom_errorbar(data = nutrientamides_emmeans, aes(x = treatment, y = response, ymin = lower.CL, ymax = upper.CL), width = 0.2, color = "black") + 
  stat_summary(fun.data = mean_se, color = "black") +
  ylab ("") +
  xlab ("Treatment")
c <- nutrientamides_graph +
  theme(legend.position = "none") + 
  scale_x_discrete(labels = c("ln+lamides" = "Low nutrients, 0.1% piperine", "hn+hamides" = "High nutrients, 2% piperine"),
                   limits = c("ln+lamides","hn+hamides"))
c

batpreference <- ggarrange(a, b, c,
                           ncol = 1, nrow = 3)
batpreference

ggsave(file="batpreference.jpg", 
       plot=batpreference,
       width=5,height=6,units="in",dpi=300)

###########################################
### For presentations, add a bat icon ####
###########################################

bat <- name_search(text = "Carollia perspicillata", options = "namebankID")[[1]] # find names
bat
bat_id_all <- name_images(uuid = bat$uid[1])  # list images
bat_id_all
bat_id <- name_images(uuid = bat$uid[1])$supertaxa[[1]]$uid  # get individual image id
bat_id
bat_pic <- image_data(bat_id, size = 256)[[1]] # get actual icon, define size. Don't run this alone
bat_pic
b_pic <- b + add_phylopic(bat_pic, alpha = 1, x = 1.5, y = 0.7, ysize = 0.4, color = "black")
b_pic
