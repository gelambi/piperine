rm(list=ls())

library(ggplot2)
library(effects)
library(glmmTMB)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(parameters)
library(viridis)
library(viridisLite)
library(rphylopic)

#################################
### Bat foraging ###
#################################

batdata <- read.csv("preference.csv")
head(batdata)
hist(batdata$percentagefoodeaten)
batdata$food_gamma <- (batdata$percentagefoodeaten/100) + 0.0001
batdata

### Filter the data per treatment. Three different two choice treatments were conducted

nutrients <- batdata %>% filter(treatment %in% c("hn", "ln"))
amides <- batdata %>% filter(treatment %in% c("ha", "la"))
nutrientamides <- batdata %>% filter(treatment %in% c("ln+lamides", "hn+hamides"))

### GLMM

# NUTRIENTS

nutrients_glmm <- glmmTMB(food_gamma ~ treatment + (1|bat), data = nutrients, family = beta_family(link="logit"))
summary(nutrients_glmm)
shapiro.test(resid(nutrients_glmm))
plot(allEffects(nutrients_glmm))
parameters(nutrients_glmm)
nutrients$prediction <- predict(nutrients_glmm, nutrients, re.form=NA,type="response")
head(nutrients)

hist(amides$percentagefoodeaten) # AMIDES
amides_glmm <- glmmTMB(food_gamma ~ treatment + (1|bat), data = amides, family = beta_family(link="logit"))
summary(amides_glmm)
shapiro.test(resid(amides_glmm))
plot(allEffects(amides_glmm))
parameters(amides_glmm)
amides$prediction <- predict(amides_glmm, amides, re.form=NA,type="response")

hist(nutrientamides$percentagefoodeaten) # NUTRIENTS AND AMIDES
nutrientamides_glmm <- glmmTMB(food_gamma ~ treatment + (1|bat), data = nutrientamides, family = beta_family(link="logit"))
summary(nutrientamides_glmm)
shapiro.test(resid(nutrientamides_glmm))
plot(allEffects(nutrientamides_glmm))
parameters(nutrientamides_glmm)
nutrientamides$prediction <- predict(nutrientamides_glmm, nutrientamides, re.form=NA,type="response")

### Graph 

nutrients_graph <- ggplot(nutrients, aes(x = treatment, y = prediction, color= treatment)) + 
  theme_classic(base_size = 13) +
  geom_boxplot(data = nutrients, aes(x = treatment, y = food_gamma, color = treatment),width = 0.1, color = "black", fill = "light grey", alpha = 0.5) + 
  geom_jitter(data = nutrients, aes(x = treatment, y = food_gamma, color = treatment), width = 0.1, size = 3, alpha = 0.7) +
  scale_color_viridis(option = "D", discrete=TRUE) +
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
  geom_boxplot(data = amides, aes(x = treatment, y = food_gamma, color = treatment),width = 0.1, color = "black", fill = "light grey", alpha = 0.5) + 
  geom_jitter(data = amides, aes(x = treatment, y = food_gamma, color = treatment), width = 0.1, size = 3, alpha = 0.7) +
  scale_color_viridis(option = "D", discrete=TRUE) +
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
  geom_boxplot(data = nutrientamides, aes(x = treatment, y = food_gamma, color = treatment),width = 0.1, color = "black", fill = "light grey", alpha = 0.5) + 
  geom_jitter(data = nutrientamides, aes(x = treatment, y = food_gamma, color = treatment), width = 0.1, size = 3, alpha = 0.7) +
  scale_color_viridis(option = "D", discrete=TRUE) +
  stat_summary(fun.data = mean_se, color = "black") +
  ylab ("") +
  xlab ("Treatment")
c <- nutrientamides_graph +
  theme(legend.position = "none") + 
  scale_x_discrete(labels = c("ln+lamides" = "Low nutrients, 0.1% piperine", "hn+hamides" = "High nutrients, 2% piperine"),
                   limits = c("ln+lamides","hn+hamides"))
c

## add a bat icon

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

batpreference <- ggarrange(a, b_pic, c,
                           ncol = 1, nrow = 3)
batpreference

ggsave(file="batpreference.jpg", 
       plot=batpreference,
       width=5,height=6,units="in",dpi=300)