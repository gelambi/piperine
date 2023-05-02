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
library(performance)

###########################################################################
### Objective 2. The effect of piperine on sugar and protein absorption ###
###########################################################################

### Absorption of proteins ###

# Read data 
setwd("~/Desktop/piperine/Nutrient absorption")
data.proteins <- read.csv("protein_absorption.csv")
data.proteins$treatment <- as.factor(data.proteins$treatment)
data.proteins$batID <- as.factor(data.proteins$batID)
data.proteins

# GLMM 

glmm1 <- glmmTMB(proteins ~ treatment + (1|date) + (1|batID), data = data.proteins)
summary(glmm1)
plot(allEffects(glmm1))
shapiro.test(resid(glmm1)) # p-value = 0.0004465 no normal. Because the residuals are not normally distributed, I have to use a different family of models

data.proteins$proteinsproportion <- data.proteins$proteins/100
data.proteins
hist(data.proteins$proteins)

glmm2 <- glmmTMB(proteinsproportion ~ treatment + (1|date) + (1|batID), data = data.proteins, beta_family(link="logit"))
summary(glmm2) # Treatment 2 0.00786 ** 
plot(allEffects(glmm2))
shapiro.test(resid(glmm2)) # p-value = 0.0004758
hist(resid(glmm2)) # they do not look too bad 
summary(allEffects(glmm2))
parameters(glmm2)
diagnose(glmm2)
r2(glmm2)
glmm2_emmeans <-emmeans(glmm2,~treatment, type="response")
glmm2_emmeans
glmm2_emmeans <- as.data.frame(glmm2_emmeans)

# Effect size 2% piperine (P = 0.008), control = 0.0159, 2% = 0.0212 
x <- 0.0212/0.0159 # 1.333333-fold
x <- 0.0212 - 0.0159
effect_size_piperine2 <- (x*100)/0.0159
effect_size_piperine2  # Percentage: the highest concentration (2%) increased protein excretion by 33.33%

## Add the model predictions
data.proteins$predictions <- predict(glmm2, data.proteins, re.form=NA,type="response")
data.proteins

protein_absorption <- ggplot(data.proteins, aes(x = treatment, y = predictions, color = treatment)) + 
  theme_classic(base_size = 11) + 
  geom_line(data = data.proteins, aes(x = treatment, y = proteinsproportion, group = batID), alpha = 0.5, color = "light grey", position = position_dodge(0.05)) +
  geom_jitter(data = data.proteins, aes(x = treatment, y = proteinsproportion), size = 2.5, position = position_dodge(0.05)) +
  scale_color_viridis(option = "D", discrete=TRUE) +
  stat_summary(fun.data = mean_se, color = "black") +
  geom_errorbar(data = glmm2_emmeans, aes(x = treatment, y = response, ymin = lower.CL, ymax = upper.CL), width = 0.2, color = "black") + 
  xlab ("Percentage of piperine in diet") +
  ylab ("Propotion of proteins excreted in fecal samples") +
  theme(legend.position = "null")

protein_absorption

ggsave(file="protein_absorption_withlines.jpg", 
       plot= protein_absorption,
       width=10,height=,units="cm",dpi=400)

### Absorption of sugars ###

# Read data 

data.sugars <- ("input_GCalignR_sugars.txt")
check_input(data = data.sugars)

# Chromatogram alignment 

peak_data_aligned <- align_chromatograms(data = data.sugars, # input data
                                         rt_col_name = "RT", # retention time variable name 
                                         rt_cutoff_low = 0, # remove peaks below 0 Minutes
                                         rt_cutoff_high = 45, # remove peaks exceeding 45 Minutes
                                         reference = NULL, # choose automatically 
                                         max_linear_shift = 0.05, # max. shift for linear corrections
                                         max_diff_peak2mean = 0.03, # max. distance of a peak to the mean across samples
                                         min_diff_peak2peak = 0.03, # min. expected distance between peaks 
                                         delete_single_peak = FALSE, # delete peaks that are present in just one sample 
                                         write_output = NULL) # add variable names to write aligned data to text files

peaksaligned <- peak_data_aligned$aligned$Area
peaksaligned
write_xlsx(peaksaligned,"output.xlsx") 

data.sugars.2 <- read.csv("output_trans.csv") # I transposed the data frame and changed the RT by letters

# The calibration curve eq. using glucose is : Area = 2736.6*Sugar concentration (ug/mL)

experiment_data <- read.csv("sugar_absorption.csv")

# Summary of calculation: Sugar concentration (ug/mL) : peak area/2736.6 = ug/mL * 1 mL (solvent) = ug * 1000 = mg / dry weight of the samples (mg) * 100 = % dry weight

data.sugars.3 <- ((((data.sugars.2[,2:18]/2736.6)/1000)/experiment_data$weight)*100)
data.sugars.3

head(experiment_data) # Now, I add some of the experiment details into the aligned peak dataframe

data.sugars.3$batID <- experiment_data$batID
data.sugars.3$treatment <- experiment_data$treatment
data.sugars.3$date <- experiment_data$date
data.sugars.3$weight <- experiment_data$weight
data.sugars.3$totalconcentration <- rowSums(data.sugars.3[,1:18])
data.sugars.3$sugarproportion <- data.sugars.3$totalconcentration/100
data.sugars.3$treatment <- as.factor(data.sugars.3$treatment)
data.sugars.3$sampleID <- experiment_data$sampleID
data.sugars.3 <- data.sugars.3[-8,] # delete an outlier
data.sugars.3

write.csv(data.sugars.3, "sugar_absorption_final.csv")
### fixed some issues in the peak table
data.sugars.3 <- read.csv("sugar_absorption_final.csv")
data.sugars.3$treatment <- as.factor(data.sugars.3$treatment)

hist(data.sugars.3$totalconcentration)
shapiro.test(data.sugars.3$totalconcentration)

glmm3 <- glmmTMB(sugarproportion ~ treatment + (1|date) + (1|batID), data = data.sugars.3)
summary(glmm3) 
plot(allEffects(glmm3))
summary(allEffects(glmm3))
shapiro.test(resid(glmm3)) # p-value = 0.0005614 no normal

glmm4 <- glmmTMB(sugarproportion ~ treatment + (1|date) + (1|batID), data = data.sugars.3, beta_family(link="logit"))
summary(glmm4) 
plot(allEffects(glmm4))
summary(allEffects(glmm4))
parameters(glmm4)
shapiro.test(resid(glmm4))
hist(resid(glmm4))
diagnose(glmm4)
glmm4_emmeans <-emmeans(glmm4,~treatment, type="response")
glmm4_emmeans <- as.data.frame(glmm4_emmeans)
glmm4_emmeans
?emmeans
r2(glmm4)
effect_size_sugars <- eff_size(glmm4_emmeans, sigma= sigma(glmm4), edf = df.residual(glmm4)) # get effect sizes 
effect_size_sugars 

## Add the model predictions

data.sugars.3$predictions <- predict(glmm4, data.sugars.3, re.form=NA,type="response")
data.sugars.3

sugar_absorption <- ggplot(data.sugars.3, aes(x = treatment, y = predictions, color = treatment)) +
  theme_classic(base_size = 11) +
  geom_line(data = data.sugars.3, aes(x = treatment, y = sugarproportion, group = batID), alpha = 0.5, color = "light grey", position = position_dodge(0.05)) +
  geom_jitter(data = data.sugars.3, aes(x = treatment, y = sugarproportion), size = 2.5, position = position_dodge(0.05)) +
  scale_color_viridis(option = "D", discrete=TRUE) +
  stat_summary(fun.data = mean_se, color = "black") +
  geom_errorbar(data = glmm4_emmeans, aes(x = treatment, y = response, ymin = lower.CL, ymax = upper.CL), width = 0.15, color = "black") +
  xlab ("Percentage of piperine in diet") +
  ylab ("Propotion of sugars excreted in fecal samples") +
  theme(legend.position = "null")

sugar_absorption

ggsave(file="sugar_absorption_withlines.jpg", 
       plot= sugar_absorption,
       width=10,height=,units="cm",dpi=400)

nutrient_absorption <- ggarrange(sugar_absorption, protein_absorption,
                        ncol = 2, nrow = 1)

nutrient_absorption

ggsave(file="nutrient_absorption_withlines.jpg", 
       plot= nutrient_absorption,
       width=20,height=10,units="cm",dpi=400)

# Multivariate analysis, NMDS: Just sugars

mmatrix <- data.sugars.3[ , 2:13] # select just the peaks
matrix <- as.matrix(mmatrix) # turn data frame into matrix

nmds_results <- metaMDS(matrix, 
                        distance = "bray",       # Specify a bray-curtis distance
                        try = 100)               # Number of iterations
nmds_results
plot(nmds_results, type = "t")

#extract NMDS scores (x and y coordinates)
data.scores <- as.data.frame(scores(nmds_results$points))
data.scores$treatment <- data.sugars.3$treatment
data.scores

#Graph

nmdsgraph <- ggplot(data.scores, aes(x = MDS1, y = MDS2, colour = treatment)) +
  geom_point(size = 5, alpha = 0.8)+
  theme_classic(base_size = 18) +
  scale_color_viridis(option = "D", discrete=TRUE, name = "Treatment (%)") +
  ylab ("NMDS2") +
  xlab ("NMDS1") + 
  stat_ellipse(level = 0.95)

nmdsgraph

ggsave(file="ind_sugar_peaks.jpg", 
       plot=nmdsgraph,
       width=8,height=5,units="in",dpi=300)

adonis <- adonis2(matrix ~ data.scores$treatment, distance = "bray", perm=9999)
adonis
write.csv(adonis, file = "adonis_sugar_peaks.csv")

######################################
### Merge proteins and sugars data ###
######################################

data.proteins <- read.csv("protein_absorption.csv")
data.proteins$treatment <- as.factor(data.proteins$treatment)
data.proteins$batID <- as.factor(data.proteins$batID)
data.proteins

data.sugars.3 <- read.csv("sugar_absorption_final.csv")
data.sugars.3$treatment <- as.factor(data.sugars.3$treatment)

bothnutrients <- merge(data.proteins, data.sugars.3, by = "sampleID")
bothnutrients <- as.data.frame(bothnutrients)
head(bothnutrients)

glmm5 <- glmmTMB(proteins ~ totalconcentration + (1|date.x) + (1|batID.x), data = bothnutrients)
summary(glmm5) 
plot(allEffects(glmm5))
hist(resid(glmm5))
shapiro.test(resid(glmm5))
summary(allEffects(glmm5))
parameters(glmm5)
diagnose(glmm5)


bothnutrients <- ggplot(bothnutrients, aes(x = totalconcentration, y = proteins, color = treatment.x)) +
  theme_classic(base_size = 13) +
  geom_jitter(width = 0.1, size = 2.5, alpha = 0.8) + 
  scale_color_viridis(option = "D", discrete=TRUE, name = "Treatment(%)") +
  xlab ("Sugars excreted in fecal samples") +
  ylab ("Proteins excreted in fecal samples")
  
bothnutrients

ggsave(file="bothnutrients.jpg", 
       plot=bothnutrients,
       width=5,height=4,units="in",dpi=300)

# With data subset 

t_0.1 <- bothnutrients %>% filter(treatment.x %in% c("0.1"))
t_0.5 <- bothnutrients %>% filter(treatment.x %in% c("0.5"))
t_1.5 <- bothnutrients %>% filter(treatment.x %in% c("1.5"))
t_2 <- bothnutrients %>% filter(treatment.x %in% c("2"))
t_0 <- bothnutrients %>% filter(treatment.x %in% c("0"))

### GLMM, is there any association between sugar and protein absorption? ###

glmm6 <- glmmTMB(proteins ~ totalconcentration + (1|date.x) + (1|batID.x), data = t_0.1)
summary(glmm6) 
hist(resid(glmm6))
shapiro.test(resid(glmm6))
parameters(glmm6)
diagnose(glmm6)
glmm6_emmeans <-emmeans(glmm6,~totalconcentration)
glmm6_emmeans
t_0.1$predictions1 <- predict(glmm6, type="response",re.form = NA, newdata = t_0.1)

bothnutrients_t_0.1 <- ggplot(t_0.1, aes(x = totalconcentration, y = proteins, color = treatment.x)) +
  theme_classic(base_size = 13) +
  geom_jitter(width = 0.1, size = 2.5, color = "#39568CFF") + 
  scale_color_viridis(option = "D", discrete=TRUE, name = "Treatment(%)") +
  geom_line(data=t_0.1, aes(x=totalconcentration,y=predictions1), linetype = 2) +
  xlab ("Excreted sugars") +
  ylab ("Excreted proteins") + 
  theme(legend.position = "null") 
bothnutrients_t_0.1

glmm7 <- glmmTMB(proteins ~ totalconcentration + (1|date.x) + (1|batID.x), data = t_0.5)
summary(glmm7) 
shapiro.test(resid(glmm7))
parameters(glmm7)
diagnose(glmm7)
glmm7_emmeans <-emmeans(glmm7,~totalconcentration)
glmm7_emmeans
t_0.5$predictions1 <- predict(glmm7, type="response",re.form = NA, newdata = t_0.5)

bothnutrients_t_0.5 <- ggplot(t_0.5, aes(x = totalconcentration, y = proteins)) +
  theme_classic(base_size = 13) +
  geom_jitter(width = 0.1, size = 2.5, color = "#238A8DFF") + 
  scale_color_viridis(option = "D", discrete=TRUE, name = "Treatment(%)") +
  geom_line(data=t_0.5, aes(x=totalconcentration,y=predictions1), linetype = 2) +
  xlab ("Excreted sugars") +
  ylab ("Excreted proteins") + 
  theme(legend.position = "null")
bothnutrients_t_0.5

glmm8 <- glmmTMB(proteins ~ totalconcentration + (1|date.x) + (1|batID.x), data = t_1.5)
summary(glmm8) 
shapiro.test(resid(glmm8))
parameters(glmm8)
diagnose(glmm8)
glmm8_emmeans <-emmeans(glmm8,~totalconcentration)
glmm8_emmeans
t_1.5$predictions1 <- predict(glmm8, type="response",re.form = NA, newdata = t_1.5)

bothnutrients_t_1.5 <- ggplot(t_1.5, aes(x = totalconcentration, y = proteins)) +
  theme_classic(base_size = 13) +
  geom_jitter(width = 0.1, size = 2.5, color = "#95D840FF" ) + 
  scale_color_viridis(option = "D", discrete=TRUE, name = "Treatment(%)") +
  geom_line(data=t_1.5, aes(x=totalconcentration,y=predictions1), linetype = 2) +
  xlab ("Excreted sugars") +
  ylab ("Excreted proteins") + 
  theme(legend.position = "null")
bothnutrients_t_1.5 

glmm9 <- glmmTMB(proteins ~ totalconcentration  + (1|date.x) + (1|batID.x), data = t_2)
summary(glmm9) 
shapiro.test(resid(glmm9))
parameters(glmm9)
diagnose(glmm9)
glmm9_emmeans <-emmeans(glmm9,~totalconcentration)
glmm9_emmeans
t_2$predictions1 <- predict(glmm9, type="response",re.form = NA, newdata = t_2)

bothnutrients_t_2 <- ggplot(t_2, aes(x = totalconcentration, y = proteins)) +
  theme_classic(base_size = 13) +
  geom_jitter(width = 0.1, size = 2.5, color = "#FDE725FF") + 
  scale_color_viridis(option = "D", discrete=TRUE, name = "Treatment(%)") +
  geom_line(data=t_2, aes(x=totalconcentration,y=predictions1), linetype = 2) +
  xlab ("Excreted sugars") +
  ylab ("Excreted proteins") + 
  theme(legend.position = "null")
bothnutrients_t_2

glmm10 <- glmmTMB(proteins ~ totalconcentration + (1|date.x) + (1|batID.x), data = t_0)
summary(glmm10) 
hist(resid(glmm10))
shapiro.test(resid(glmm10))
parameters(glmm10)
diagnose(glmm10)
glmm10_emmeans <-emmeans(glmm10,~totalconcentration)
glmm10_emmeans
t_0$predictions1 <- predict(glmm10, type="response",re.form = NA, newdata = t_0)

bothnutrients_t_0 <- ggplot(t_0, aes(x = totalconcentration, y = proteins, color = treatment.x)) +
  theme_classic(base_size = 13) +
  geom_jitter(width = 0.1, size = 2.5) + 
  scale_color_viridis(option = "D", discrete=TRUE, name = "Treatment(%)") +
  geom_line(data=t_0, aes(x=totalconcentration,y=predictions1), linetype = 2) +
  xlab ("Excreted sugars") +
  ylab ("Excreted proteins") + 
  theme(legend.position = "null") 
bothnutrients_t_0

ind_molecules_2 <- ggarrange(bothnutrients_t_0,
                             bothnutrients_t_0.1,
                             bothnutrients_t_0.5,
                             bothnutrients_t_1.5,
                             bothnutrients_t_2,
                             ncol = 3, nrow = 2)
ind_molecules_2

ggsave(file="ind_molecules_2.jpg", 
       plot= ind_molecules_2,
       width=15,height=12,units="cm",dpi=400)

