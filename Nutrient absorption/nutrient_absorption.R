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
library(car)
library(gplots)
library(reshape2)
library(tidyr)

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
data.proteins <- data.proteins %>% filter(!batID %in% c("1"))

# GLMM 

glmm1 <- glmmTMB(proteins ~ treatment + (1|date) + (1|batID), data = data.proteins)
summary(glmm1)
plot(allEffects(glmm1))
shapiro.test(resid(glmm1)) # p-value = 0.0004465 no normal. Because the residuals are not normally distributed, I have to use a different family of models

data.proteins$proteinsproportion <- data.proteins$proteins/100
data.proteins
hist(data.proteins$proteins)

glmm2 <- glmmTMB(proteinsproportion ~ treatment + (1|date) + (1|batID), data = data.proteins, beta_family(link="logit"))
summary(glmm2) # Treatment 2 0.00955 ** 
Anova(glmm2)
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

# Effect size 2% piperine (P = 0.008), control = 0.0161, 2% = 0.0213 

x <- 0.0213/0.0161
x  # 1.3-fold
x <- 0.0213 - 0.0161
x # 0.0052
effect_size_piperine <- (x*100)/0.0151
effect_size_piperine  # Percentage: the highest concentration (2%) increased protein excretion by 34%

## Add the model predictions
data.proteins$predictions <- predict(glmm2, data.proteins, re.form=NA,type="response")
data.proteins

protein_absorption <- ggplot(data.proteins, aes(x = treatment, y = predictions, color = batID)) + 
  theme_classic(base_size = 11) + 
  #geom_line(data = data.proteins, aes(x = treatment, y = proteinsproportion, group = batID), alpha = 0.3, color = "light grey", position = position_dodge(0.05)) +
  geom_violin(data = data.proteins, aes(x = treatment, y = proteinsproportion), fill = "light gray", color = "white", alpha = 0.4) +
  geom_jitter(data = data.proteins, aes(x = treatment, y = proteinsproportion), size = 2.5, position = position_dodge(0.05)) +
  scale_color_viridis(option = "D", discrete=TRUE) +
  stat_summary(fun.data = mean_se, color = "black") +
  geom_errorbar(data = glmm2_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "black") + 
  xlab ("Percentage of piperine in diet") +
  ylab ("Proportion of total proteins excreted") +
  theme(legend.position = "null") +
  ggtitle("(A)") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_text(aes(x = "2", y = 0.06, label = "**"), size = 5, color = "black")

protein_absorption

ggsave(file="protein_absorption_withlines.jpg", 
       plot= protein_absorption,
       width=10,height=,units="cm",dpi=400)

### as numeric 
data.proteins$treatment  <- as.numeric(data.proteins$treatment)
glmm2 <- glmmTMB(proteinsproportion ~ treatment + (1|date) + (1|batID), data = data.proteins)
summary(glmm2) # Treatment 2 0.00955 ** 
plot(allEffects(glmm2))
shapiro.test(resid(glmm2))

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
setwd("/Users/marianagelambi/Desktop/piperine/Nutrient absorption")
data.sugars.3 <- read.csv("sugar_absorption_final.csv")
data.sugars.3$treatment <- as.factor(data.sugars.3$treatment)
data.sugars.3$batID <- as.factor(data.sugars.3$batID)
data.sugars.3 <- data.sugars.3 %>% filter(!batID %in% c("1"))
str(data.sugars.3)

data.sugars.3
justpeaks <- data.sugars.3[, 2:13]

summary_indsugars <- data.frame(
  Variable = names(justpeaks),
  Mean = sapply(justpeaks, function(x) mean(x, na.rm = TRUE)),
  SD = sapply(justpeaks, function(x) sd(x, na.rm = TRUE)),
  Freq_Detection = sapply(justpeaks, function(x) sum(x != 0, na.rm = TRUE)),
  Percentage_Detection = sapply(justpeaks, function(x) (sum(x != 0, na.rm = TRUE) / length(x)) * 100)
)
write.csv(summary_indsugars, "summary_indsugars.csv")

### create a representative chromatogram: 

draw_chromatogram <- data.frame(
  mean_RT = c(1.39, 1.825, 3.503, 3.803, 4.253),
  peak_11 = c(3.29, 0.1, 0.46, 0.63, 0.40), # actually the mean abundance, but it works
  Width = c(0.05, 0.05, 0.05, 0.05, 0.05) # it is not informative
)

# Create a data frame to store the Gaussian curves
gaussian_curves <- data.frame()
n_points <- 1000  # Number of points for each curve

for (i in 1:nrow(draw_chromatogram)) {
  curve_data <- data.frame(
    x = seq(draw_chromatogram$mean_RT[i] - 3 * draw_chromatogram$Width[i],
            draw_chromatogram$mean_RT[i] + 3 * draw_chromatogram$Width[i], length.out = n_points),
    y = dnorm(
      seq(draw_chromatogram$mean_RT[i] - 3 * draw_chromatogram$Width[i],
          draw_chromatogram$mean_RT[i] + 3 * draw_chromatogram$Width[i], length.out = n_points),
      mean = draw_chromatogram$mean_RT[i],
      sd = draw_chromatogram$Width[i]
    ) * draw_chromatogram$peak_11[i]
  )
  gaussian_curves <- rbind(gaussian_curves, curve_data)
}

# Sum all the Gaussian components to create a smooth line
smooth_line <- gaussian_curves %>%
  group_by(x) %>%
  summarise(y = sum(y))

# Create the chromatogram plot
representative_chromatogram <- ggplot() +
  theme_pubr() + 
  geom_line(data = smooth_line, aes(x = x, y = y), linewidth = 0.9, color = "darkmagenta") +
  labs(x = "Retention Time (min)", y = "Intensity") + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  geom_text(aes(x = 1.39, y = 27, label = "A"), size = 5, color = "black") + 
  geom_text(aes(x = 1.825, y = 2, label = "B"), size = 5, color = "black") +
  geom_text(aes(x = 3.503, y = 5, label = "F"), size = 5, color = "black") + 
  geom_text(aes(x = 3.803, y = 6.2, label = "G"), size = 5, color = "black") +
  geom_text(aes(x = 4.253, y = 4.5, label = "K"), size = 5, color = "black")

representative_chromatogram
ggsave(file="representative_chromatogram.jpg", 
       plot= representative_chromatogram,
       width=10,height=10,units="cm",dpi=400)
#########

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
r2(glmm4)

## Add the model predictions

data.sugars.3$predictions <- predict(glmm4, data.sugars.3, re.form=NA,type="response")

sugar_absorption <- ggplot(data.sugars.3, aes(x = treatment, y = predictions, color = batID)) +
  theme_classic(base_size = 11) +
  #geom_line(data = data.sugars.3, aes(x = treatment, y = sugarproportion, group = batID), alpha = 0.3, color = "light grey", position = position_dodge(0.05)) +
  geom_violin(data = data.sugars.3, aes(x = treatment, y = sugarproportion),  fill = "light gray", color = "white", alpha = 0.4) +
   geom_jitter(data = data.sugars.3, aes(x = treatment, y = sugarproportion), size = 2.5, position = position_dodge(0.05)) +
  scale_color_viridis(option = "D", discrete=TRUE) +
  stat_summary(fun.data = mean_se, color = "black") +
  geom_errorbar(data = glmm4_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15, color = "black") +
  xlab ("Percentage of piperine in diet") +
  ylab ("Proportion of total sugars excreted") +
  theme(legend.position = "none") + 
  ggtitle("(B)") +
  theme(plot.title = element_text(hjust = 0.5))

sugar_absorption

ggsave(file="sugar_absorption_withlines.jpg", 
       plot= sugar_absorption,
       width=10,height=,units="cm",dpi=400)

nutrient_absorption <- ggarrange(protein_absorption,
                                 sugar_absorption,
                                 ncol = 2, nrow = 1)

nutrient_absorption

ggsave(file="nutrient_absorption_withlines.jpg", 
       plot= nutrient_absorption,
       width=20,height=10,units="cm",dpi=400)

### Individual sugars
data.sugars.3$A
data.sugars.3$A_beta <- (data.sugars.3$A)/100
glmm5 <- glmmTMB(A_beta ~ treatment + (1|date) + (1|batID), data = data.sugars.3, beta_family(link="logit"))
summary(glmm5) 
Anova(glmm5)
plot(allEffects(glmm5))
summary(allEffects(glmm5))
parameters(glmm5)
shapiro.test(resid(glmm5))
diagnose(glmm5)
r2(glmm5)
glmm5_emmeans <-emmeans(glmm5,~treatment, type="response")
glmm5_emmeans
glmm5_emmeans <- as.data.frame(glmm5_emmeans)

data.sugars.3$predictions <- predict(glmm5, data.sugars.3, re.form=NA,type="response")

sugar_absorption_A <- ggplot(data.sugars.3, aes(x = treatment, y = predictions, color = batID)) +
  theme_classic(base_size = 10) +
  #geom_violin(data = data.sugars.3, aes(x = treatment, y = A_beta),  fill = "light gray", color = "white", alpha = 0.4) +
  geom_jitter(data = data.sugars.3, aes(x = treatment, y = A_beta), alpha = 0.5, size = 2, position = position_dodge(0.05)) +
  scale_color_viridis(option = "D", discrete=TRUE) +
  stat_summary(fun.data = mean_se, color = "black") +
  geom_errorbar(data = glmm5_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15, color = "black") +
  xlab (" ") +
  ylab ("A") +
  theme(legend.position = "none") + 
  geom_text(aes(x = "2", y = 0.075, label = "*"), size = 8, color = "black")
sugar_absorption_A

# Effect size 2% piperine (P = 0.008), control = 0.0161, 2% = 0.0213 

glmm5_emmeans
x <- 0.0319/0.0290
x  # 1.1-fold
x <- 0.0319-0.0290
x # 0.0029
effect_size_piperine <- (x*100)/0.0319
effect_size_piperine 

data.sugars.3$B
data.sugars.3$B_beta <- ((data.sugars.3$B)/100) + 0.001
glmm6 <- glmmTMB(B_beta ~ treatment + (1|date) + (1|batID), data = data.sugars.3, beta_family(link="logit"))
summary(glmm6) 
Anova(glmm6)
plot(allEffects(glmm6))
summary(allEffects(glmm6))
parameters(glmm6)
shapiro.test(resid(glmm6))
diagnose(glmm6)
r2(glmm6)
glmm6_emmeans <-emmeans(glmm6,~treatment, type="response")
glmm6_emmeans
glmm6_emmeans <- as.data.frame(glmm6_emmeans)

data.sugars.3$predictions <- predict(glmm6, data.sugars.3, re.form=NA,type="response")

sugar_absorption_B <- ggplot(data.sugars.3, aes(x = treatment, y = predictions, color = batID)) +
  theme_classic(base_size = 10) +
  #geom_violin(data = data.sugars.3, aes(x = treatment, y = B_beta),  fill = "light gray", color = "white", alpha = 0.4) +
  geom_jitter(data = data.sugars.3, aes(x = treatment, y = B_beta), alpha = 0.5, size = 2, position = position_dodge(0.05)) +
  scale_color_viridis(option = "D", discrete=TRUE) +
  stat_summary(fun.data = mean_se, color = "black") +
  geom_errorbar(data = glmm6_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15, color = "black") +
  xlab (" ") +
  ylab ("B") +
  theme(legend.position = "none")
sugar_absorption_B

data.sugars.3$F
data.sugars.3$F_beta <- ((data.sugars.3$F)/100) + 0.001
glmm7 <- glmmTMB(F_beta ~ treatment + (1|date) + (1|batID), data = data.sugars.3, beta_family(link="logit"))
summary(glmm7) 
Anova(glmm7)
plot(allEffects(glmm7))
summary(allEffects(glmm7))
parameters(glmm7)
shapiro.test(resid(glmm7))
diagnose(glmm7)
r2(glmm7)
glmm7_emmeans <-emmeans(glmm7,~treatment, type="response")
glmm7_emmeans
glmm7_emmeans <- as.data.frame(glmm7_emmeans)

data.sugars.3$predictions <- predict(glmm7, data.sugars.3, re.form=NA,type="response")

sugar_absorption_F <- ggplot(data.sugars.3, aes(x = treatment, y = predictions, color = batID)) +
  theme_classic(base_size = 10) +
  #geom_violin(data = data.sugars.3, aes(x = treatment, y = F_beta),  fill = "light gray", color = "white", alpha = 0.4) +
  geom_jitter(data = data.sugars.3, aes(x = treatment, y = F_beta), alpha = 0.5, size = 2, position = position_dodge(0.05)) +
  scale_color_viridis(option = "D", discrete=TRUE) +
  stat_summary(fun.data = mean_se, color = "black") +
  geom_errorbar(data = glmm7_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15, color = "black") +
  xlab ("") +
  ylab ("F") +
  theme(legend.position = "none")
sugar_absorption_F

data.sugars.3$G
data.sugars.3$G_beta <- ((data.sugars.3$G)/100) + 0.001
glmm8 <- glmmTMB(G_beta ~ treatment + (1|date) + (1|batID), data = data.sugars.3, beta_family(link="logit"))
summary(glmm8) 
Anova(glmm8)
plot(allEffects(glmm8))
summary(allEffects(glmm8))
parameters(glmm8)
shapiro.test(resid(glmm8))
diagnose(glmm8)
r2(glmm8)
glmm8_emmeans <-emmeans(glmm8,~treatment, type="response")
glmm8_emmeans
glmm8_emmeans <- as.data.frame(glmm8_emmeans)

glmm8_emmeans
x <- 0.00796-0.00605
x # 0.00191
effect_size_piperine <- (x*100)/0.00796
effect_size_piperine 

data.sugars.3$predictions <- predict(glmm8, data.sugars.3, re.form=NA,type="response")

sugar_absorption_G <- ggplot(data.sugars.3, aes(x = treatment, y = predictions, color = batID)) +
  theme_classic(base_size = 10) +
  #geom_violin(data = data.sugars.3, aes(x = treatment, y = G_beta),  fill = "light gray", color = "white", alpha = 0.4) +
  geom_jitter(data = data.sugars.3, aes(x = treatment, y = G_beta), alpha = 0.5, size = 2, position = position_dodge(0.05)) +
  scale_color_viridis(option = "D", discrete=TRUE) +
  stat_summary(fun.data = mean_se, color = "black") +
  geom_errorbar(data = glmm8_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15, color = "black") +
  xlab ("Piperine (%)") +
  ylab ("G") +
  theme(legend.position = "none") + 
  geom_text(aes(x = "1.5", y = 0.02, label = "*"), size = 8, color = "black")
sugar_absorption_G

data.sugars.3$K
data.sugars.3$K_beta <- ((data.sugars.3$K)/100) + 0.001
glmm9 <- glmmTMB(K_beta ~ treatment + (1|date) + (1|batID), data = data.sugars.3, beta_family(link="logit"))
summary(glmm9) 
Anova(glmm9)
plot(allEffects(glmm9))
summary(allEffects(glmm9))
parameters(glmm9)
shapiro.test(resid(glmm9))
diagnose(glmm9)
r2(glmm9)
glmm9_emmeans <-emmeans(glmm9,~treatment, type="response")
glmm9_emmeans
glmm9_emmeans <- as.data.frame(glmm9_emmeans)

data.sugars.3$predictions <- predict(glmm9, data.sugars.3, re.form=NA,type="response")

sugar_absorption_K <- ggplot(data.sugars.3, aes(x = treatment, y = predictions, color = batID)) +
  theme_classic(base_size = 10) +
  #geom_violin(data = data.sugars.3, aes(x = treatment, y = K_beta),  fill = "light gray", color = "white", alpha = 0.4) +
  geom_jitter(data = data.sugars.3, aes(x = treatment, y = K_beta), alpha = 0.5, size = 2, position = position_dodge(0.05)) +
  scale_color_viridis(option = "D", discrete=TRUE) +
  stat_summary(fun.data = mean_se, color = "black") +
  geom_errorbar(data = glmm9_emmeans, aes(x = treatment, y = response, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15, color = "black") +
  xlab ("Piperine (%)") +
  ylab ("K") +
  theme(legend.position = "none") + 
  geom_text(aes(x = "2", y = 0.011, label = "*"), size = 8, color = "black")
sugar_absorption_K

glmm9_emmeans
x <- 0.004166173/0.002882884
x  # 1.4-fold
x <- 0.004166173-0.002882884
x # 0.001283289
effect_size_piperine <- (x*100)/0.004166173
effect_size_piperine 

individual_sugars <- ggarrange(sugar_absorption_A,
                               sugar_absorption_B,
                               sugar_absorption_F,
                               sugar_absorption_G,
                               sugar_absorption_K,
                               align = "hv",
                               ncol = 2, nrow = 3)
individual_sugars

ggsave(file="individual_sugars.jpg", 
       plot=individual_sugars,
       width=5,height=5,units="in",dpi=300)


### before NMDS, evaluate correlation between sugar peaks
sugar_peaks <- data.sugars.3[, c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L")]
# Calculate the correlation matrix
cor_matrix <- cor(sugar_peaks)
p.mat <- cor.mtest(cor_matrix)
p.mat <- p.mat$p
# Define significance levels
significance_levels <- c(0.05, 0.01, 0.001)
significance_labels <- c("*", "**", "***")
# Function to add significance labels
add_significance_labels <- function(p.mat, significance_levels, significance_labels) {
  sig_matrix <- matrix("", nrow = ncol(p.mat), ncol = ncol(p.mat))
  for (i in 1:length(significance_levels)) {
    sig_matrix[p.mat <= significance_levels[i]] <- significance_labels[i]
  }
  return(sig_matrix)
}
sig_labels_matrix <- add_significance_labels(p.mat, significance_levels, significance_labels)

# Create a data frame for the heatmap
cor_df <- as.data.frame(cor_matrix)
cor_df <- cor_df %>% 
  rownames_to_column(var = "Peaks1") %>%
  melt(id.vars = "Peaks1", variable.name = "Peaks2", value.name = "Correlation")

# Create a data frame for significance labels
sig_labels_df <- as.data.frame(sig_labels_matrix)
sig_labels_df <- sig_labels_df %>% 
  rownames_to_column(var = "Peaks1") %>%
  melt(id.vars = "Peaks1", variable.name = "Peaks2", value.name = "Label")
# Complete the data frame to match the dimensions of the heatmap
complete_df <- complete(cor_df, Peaks1, Peaks2)
complete_df$Label <- sig_labels_df$Label

correlation_sugars <- ggplot(complete_df, aes(x = Peaks1, y = Peaks2, fill = Correlation)) +
  geom_tile() +
  labs(x = " ",
       y = " ") +
  scale_fill_gradient2(high =  "#400154BF", mid = "white", low = "#D8E219FF",
                       midpoint = 0, limits = c(-1, 1), name = "Correlation\nCoefficient") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = Label), color = "black", size = 4) +
  coord_fixed(ratio = 1)
correlation_sugars 

ggsave(file="correlation_sugars.jpg", 
       plot=correlation_sugars,
       width=5,height=5,units="in",dpi=300)

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
data.scores$batID <- data.sugars.3$batID
data.scores

#Graph
color_palette <- viridis(20, option = "A")
selected_colors <- color_palette[c(3, 7, 11, 14, 17)]

library(MoMAColors)

nmdsgraph <- ggplot(data.scores, aes(x = MDS1, y = MDS2, color = treatment)) +
  geom_point(size = 4, alpha = 0.8)+
  theme_classic(base_size = 18) +
  scale_color_moma_d("Kippenberger", name = "Treatment (%)") +
  ylab ("NMDS2") +
  xlab ("NMDS1") + 
  stat_ellipse(level = 0.95)
nmdsgraph

ggsave(file="ind_sugar_peaks.jpg", 
       plot=nmdsgraph,
       width=6,height=4,units="in",dpi=300)

# PERMANOVA
adonis <- adonis2(matrix ~ data.scores$treatment, distance = "bray", perm=9999)
adonis
write.csv(adonis, file = "adonis_sugar_peaks.csv")
# BETADISP
dist_matrix <- vegdist(matrix, method = "bray")
groups <- data.sugars.3$treatment
dispersal <- betadisper(dist_matrix, groups, type = c("centroid"))
anova(dispersal)
