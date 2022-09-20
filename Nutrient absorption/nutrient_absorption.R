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

### Absorption of proteins ###

# Read data 

data.proteins <- read.csv("protein_absorption.csv")
data.proteins$treatment <- as.factor(data.proteins$treatment)
data.proteins$batID <- as.factor(data.proteins$batID)
data.proteins

# GLMM 

glmm1 <- glmmTMB(proteins ~ treatment + (1|date) + (1|batID), data = data.proteins)
summary(glmm1)
shapiro.test(resid(glmm1)) # p-value = 0.0004465 no normal
plot(allEffects(glmm1))

# Because the residuals are not normally distributed, I have to use a different family of models

glmm2 <- glmmTMB(proteins ~ treatment + (1|date) + (1|batID), data = data.proteins, family= "Gamma") # I cannot run this model because proteins is in %, so I would create a new column of the same values but as proportion

data.proteins$proteinsproportion <- data.proteins$proteins/100
data.proteins

hist(data.proteins$proteins)

glmm2 <- glmmTMB(proteinsproportion ~ treatment + (1|date) + (1|batID), data = data.proteins, family=Gamma(link = log))
summary(glmm2) 
plot(allEffects(glmm2))
shapiro.test(resid(glmm2))
summary(allEffects(glmm2)) ## Treatment 2% is significantly different from the control (intercept) p = 0.00786.
parameters(glmm2)

glmm2 <- glmmTMB(proteinsproportion ~ treatment + (1|date) + (1|batID), data = data.proteins, beta_family(link="logit"))
summary(glmm2) 
plot(allEffects(glmm2))
shapiro.test(resid(glmm2))
summary(allEffects(glmm2)) ## Treatment 2% is significantly different from the control (intercept) p = 0.00786.
parameters(glmm2)

## Add the model predictions

data.proteins$predictions <- predict(glmm2, data.proteins, re.form=NA,type="response")
data.proteins

protein_absorption <- ggplot(data.proteins, aes(x = treatment, y = predictions, color = treatment)) +
  theme_classic(base_size = 11) +
  geom_boxplot(data = data.proteins, aes(x = treatment, y = proteinsproportion), alpha = 0.8, color = "grey", fill = "light grey", width = 0.3) +
  geom_jitter(data = data.proteins, aes(x = treatment, y = proteinsproportion, color = treatment), width = 0.1, size = 2.5, alpha = 0.8) +
  scale_color_viridis(option = "D", discrete=TRUE) +
  stat_summary(fun.data = mean_se, color = "black") +
  xlab ("Percentage of piperine in diet") +
  ylab ("Propotion of proteins excreted in fecal samples") +
  theme(legend.position = "null")

protein_absorption

ggsave(file="protein_absorption.jpg", 
       plot= protein_absorption,
       width=10,height=,units="cm",dpi=400)


##################

### Absorption of sugars ###

# Read data 

data.sugars <- ("input_GCalignR_sugars.txt")
check_input(data = data.sugars)

# Chromatogram aligment 

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

data.sugars.2 <- read.csv("output_trans.csv") # I transposed the dataframe and changed the RT by letters

# The calibration curve eq. using glucose is : Area = 2736.6*Sugar concentration (ug/mL)

experiment_data <- read.csv("sugar_absorption.csv")

# Summary of calculation: Sugar concentration (ug/mL) : peak area/2736.6 = ug/mL * 1 mL (solvent) = ug * 1000 = mg / dry weight of the samples (mg) * 100 = % dry weight

data.sugars.3 <- ((((data.sugars.2[,2:18]/2736.6)/1000)/experiment_data$weight)*100)
data.sugars.3

head(experiment_data) # Now, I add some of the experiment details into the aligned peak dataframe

data.sugars.3$batID <- experiment_data$batTD
data.sugars.3$treatment <- experiment_data$treatment
data.sugars.3$date <- experiment_data$date
data.sugars.3$weight <- experiment_data$weight
data.sugars.3$totalconcentration <- rowSums(data.sugars.3[,1:18])
data.sugars.3$sugarproportion <- data.sugars.3$totalconcentration/100
data.sugars.3$treatment <- as.factor(data.sugars.3$treatment)
data.sugars.3 <- data.sugars.3[-8,] # delete an outlier
data.sugars.3

hist(data.sugars.3$totalconcentration)

glmm3 <- glmmTMB(sugarproportion ~ treatment + (1|date) + (1|batID), data = data.sugars.3)
summary(glmm3) 
plot(allEffects(glmm3))
summary(allEffects(glmm3))
shapiro.test(resid(glmm3)) # p-value = 0.0005614 no normal

glmm4 <- glmmTMB(sugarproportion ~ treatment + (1|date) + (1|batID), data = data.sugars.3, family=Gamma(link = log))
summary(glmm4) 
plot(allEffects(glmm4))
summary(allEffects(glmm4))
parameters(glmm4)

## Add the model predictions

data.sugars.3$predictions <- predict(glmm4, data.sugars.3, re.form=NA,type="response")
data.sugars.3

sugar_absorption <- ggplot(data.sugars.3, aes(x = treatment, y = predictions, color = treatment)) +
  theme_classic(base_size = 11) +
  geom_boxplot(data = data.sugars.3, aes(x = treatment, y = sugarproportion), alpha = 0.8, color = "grey", fill = "light grey", width = 0.3) +
  geom_jitter(data = data.sugars.3, aes(x = treatment, y = sugarproportion, color = treatment), width = 0.1, size = 2.5, alpha = 0.8) +
  scale_color_viridis(option = "D", discrete=TRUE) +
  stat_summary(fun.data = mean_se, color = "black") +
  xlab ("Percentage of piperine in diet") +
  ylab ("Propotion of sugars excreted in fecal samples") +
  theme(legend.position = "null")

sugar_absorption

ggsave(file="sugar_absorption.jpg", 
       plot= sugar_absorption,
       width=10,height=,units="cm",dpi=400)

nutrient_absorption <- ggarrange(sugar_absorption, protein_absorption,
                        ncol = 2, nrow = 1)

nutrient_absorption

ggsave(file="nutrient_absorption.jpg", 
       plot= nutrient_absorption,
       width=20,height=10,units="cm",dpi=400)

# Multivariate analysis, NMDS: Just sugars

mmatrix <- data.sugars.3[ , 1:18] # select just the peaks
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
  xlab ("NMDS1") 

nmdsgraph

ggsave(file="ind_sugar_peaks.jpg", 
       plot=nmdsgraph,
       width=8,height=5,units="in",dpi=300)

adonis <- adonis2(matrix ~ data.scores$treatment, distance = "bray", perm=9999)
adonis
write.csv(adonis, file = "adonis_sugar_peaks.csv")

# Now, I will try to see if the treatment affected unique peaks 

glmm_I <- glmmTMB(I ~ treatment + (1|date) + (1|batID), data = data.sugars.3)
summary(glmm_I) 
plot(allEffects(glmm_I))
summary(allEffects(glmm_I))

glmm_A <- glmmTMB(A ~ treatment + (1|date) + (1|batID), data = data.sugars.3)
summary(glmm_A) 
plot(allEffects(glmm_A))
summary(allEffects(glmm_A))

glmm_K <- glmmTMB(K ~ treatment + (1|date) + (1|batID), data = data.sugars.3)
summary(glmm_K) 
plot(allEffects(glmm_K))
summary(allEffects(glmm_K))
