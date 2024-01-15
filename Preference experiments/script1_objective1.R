#################################################################################################
### Objective 1. The relative role of nutrients and defensive metabolites in bat preference ###
#################################################################################################

rm(list=ls()) 

library(tidyverse)
library(ggplot2)
library(effects)
library(ggplot2)
library(glmmTMB)
library(viridis)
library(viridisLite)
library(RColorBrewer)
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
library(sjPlot)
library(performance)

### Objectives:
### (1.1) Do bats distinguish between a high- and low-nutrient diet?
### (1.2) Do bats avoid fruit secondary metabolites?
### (1.3) Does a high-nutrient diet mitigate the deterrent effects of high concentrations of secondary metabolites, as postulated by the 'nutrient-toxin titration hypothesis’?
### (1.4) Do bats adjust their overall nutrient or secondary metabolite intake based on the food options available?

#################
### READ DATA ###
#################

setwd("~/Desktop/piperine/Preference experiments")
batdata <- read.csv("preference.csv")
batdata$bat <- as.factor(batdata$bat)
hist(batdata$percentagefoodeaten)
batdata

#####################
### Paired t-test ###
#####################

# nutrients
nutrients <- batdata %>%
  filter(treatment %in% c("hn", "ln"))
head(nutrients)
nutrients_ttest <- nutrients[, c("bat", "treatment", "foodconsumed")]

nutrients_hn <- nutrients_ttest %>%
  filter(treatment %in% c("hn"))
nutrients_ln <- nutrients_ttest %>%
  filter(treatment %in% c("ln"))

nutrients_ttest_total <- left_join(nutrients_hn, nutrients_ln, by = "bat")
# Replace values less than 0 with 0
nutrients_ttest_total$foodconsumed.x <- ifelse(nutrients_ttest_total$foodconsumed.x < 0, 0, nutrients_ttest_total$foodconsumed.x)
nutrients_ttest_total$foodconsumed.y <- ifelse(nutrients_ttest_total$foodconsumed.y < 0, 0, nutrients_ttest_total$foodconsumed.y)

result_nutrients <- t.test(nutrients_ttest_total$foodconsumed.x, nutrients_ttest_total$foodconsumed.y, mu = 0, paired = TRUE)
result_nutrients
head(nutrients_ttest_total)
nutrients_ttest_total$diff <- nutrients_ttest_total$foodconsumed.x - nutrients_ttest_total$foodconsumed.y
mean(nutrients_ttest_total$diff) # this is the mean difference in the t-test
2.9063169 - 0.4193498

nutrients_ttest_total$difference <- (nutrients_ttest_total$foodconsumed.x - nutrients_ttest_total$foodconsumed.y)

# Create a horizontal line plot with colors for positive and negative values
plot_nutrients <- ggplot(nutrients_ttest_total, aes(y = difference)) +
  theme_test(base_size = 10) + 
  ggtitle("(A) Trial 1.1")+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "bottom") + 
  geom_violin(aes(x = 0, y = difference), fill = "light gray", color = "white", alpha = 0.4) + 
  geom_point(aes(x = 0, color = factor(difference > 0)), size = 5, alpha = 0.7) +
  annotate("segment", x = -Inf, xend = Inf, y = 0, yend = 0, color = "black", linetype = "dashed") + 
  xlab (" ") + 
  ylab ("Differences in the amount of diet consumed (g)\nHigh nutrients − Low nutrients") + 
  scale_x_continuous(breaks = NULL) +  # Remove x-axis tick marks
  guides(color = guide_legend(title = NULL)) +  # Remove legend title
  scale_color_manual(values = c("TRUE" = "#440154FF", "FALSE" = "darkorange"),
                     labels = c("Low nutrient\npreferred", "High nutrient\npreferred"),
                     name = "Preference") + 
  scale_y_continuous(breaks = c(-5, 0, 5),
                     limits = c(-5, 5))  # Set the limits of the Y-axis
plot_nutrients

### “Geometric Framework for Nutrition”

### calculate the amount of nutrients (protein and carbohydrates).
nutrients_ttest_total$nutrients_dietX <- ((nutrients_ttest_total$foodconsumed.x)*179/5000)*1000
nutrients_ttest_total$nutrients_dietY <- ((nutrients_ttest_total$foodconsumed.y)*65/5000)*1000
nutrients_ttest_total$nutrient_total <- nutrients_ttest_total$nutrients_dietX + nutrients_ttest_total$nutrients_dietY

t1.1_nutrients <- c(nutrients_ttest_total$nutrient_total)

# piperine 
piperine <- batdata %>%
  filter(treatment %in% c("ha", "la"))
head(piperine)
piperine_ttest <- piperine[, c("bat", "treatment", "foodconsumed")]

piperine_hn <- piperine_ttest %>%
  filter(treatment %in% c("ha"))
piperine_ln <- piperine_ttest %>%
  filter(treatment %in% c("la"))

piperine_ttest_total <- left_join(piperine_hn, piperine_ln, by = "bat")
# Replace values less than 0 with 0

piperine_ttest_total$foodconsumed.x <- ifelse(piperine_ttest_total$foodconsumed.x < 0, 0, piperine_ttest_total$foodconsumed.x)
piperine_ttest_total$foodconsumed.y <- ifelse(piperine_ttest_total$foodconsumed.y < 0, 0, piperine_ttest_total$foodconsumed.y)

result_piperine <- t.test(piperine_ttest_total$foodconsumed.x, piperine_ttest_total$foodconsumed.y, mu = 0, paired = TRUE)
result_piperine

piperine_ttest_total$difference <- (piperine_ttest_total$foodconsumed.x - piperine_ttest_total$foodconsumed.y)

# Create a horizontal line plot with colors for positive and negative values
plot_piperine <- ggplot(piperine_ttest_total, aes(y = difference)) +
  theme_test(base_size = 10) + 
  ggtitle("(B) Trial 1.2") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "bottom") +
  geom_violin(aes(x = 0, y = difference), fill = "light gray", color = "white", alpha = 0.4) + 
  geom_point(aes(x = 0, color = factor(difference > 0)), size = 5, alpha = 0.7) +
  annotate("segment", x = -Inf, xend = Inf, y = 0, yend = 0, color = "black", linetype = "dashed") + 
  xlab (" ") + 
  ylab ("2% piperine − 0.1% piperine") + 
  scale_x_continuous(breaks = NULL) +  # Remove x-axis tick marks
  guides(color = guide_legend(title = NULL)) +  # Remove legend title
  scale_color_manual(values = c("TRUE" = "#440154FF", "FALSE" = "darkorange"),
                     labels = c("0.1% piperine\npreferred", "2% piperine\npreferred"),
                     name = "Preference") + 
  scale_y_continuous(breaks = c(-5, 0, 5),
                     limits = c(-5, 5))  # Set the limits of the Y-axis
plot_piperine

### “Geometric Framework for Nutrition”

head(piperine_ttest_total)
piperine_ttest_total$nutrients_dietX <- ((piperine_ttest_total$foodconsumed.x)*179/5000)*1000
piperine_ttest_total$nutrients_dietY <- ((piperine_ttest_total$foodconsumed.y)*179/5000)*1000
piperine_ttest_total$nutrient_total <- piperine_ttest_total$nutrients_dietX + piperine_ttest_total$nutrients_dietY
t1.2_nutrients <- c(piperine_ttest_total$nutrient_total)

piperine_ttest_total$piperine_dietX <- ((piperine_ttest_total$foodconsumed.x)*14.5/5000)*1000
piperine_ttest_total$piperine_dietY <- ((piperine_ttest_total$foodconsumed.y)*0.7/5000)*1000
piperine_ttest_total$piperine_total <- piperine_ttest_total$piperine_dietX + piperine_ttest_total$piperine_dietY
t1.2_piperine <- c(piperine_ttest_total$piperine_total)

# nutrients + piperine 
nutrients_piperine <- batdata %>%
  filter(treatment %in% c("ln+lamides", "hn+hamides"))
head(nutrients_piperine)
nutrients_piperine_ttest <- nutrients_piperine[, c("bat", "treatment", "foodconsumed")]

nutrients_piperine_ln <- nutrients_piperine_ttest %>%
  filter(treatment %in% c("ln+lamides"))
nutrients_piperine_hn <- nutrients_piperine_ttest %>%
  filter(treatment %in% c("hn+hamides"))

nutrients_piperine_ttest_total <- left_join(nutrients_piperine_hn, nutrients_piperine_ln, by = "bat")
# Replace values less than 0 with 0
nutrients_piperine_ttest_total$foodconsumed.x <- ifelse(nutrients_piperine_ttest_total$foodconsumed.x < 0, 0, nutrients_piperine_ttest_total$foodconsumed.x)
nutrients_piperine_ttest_total$foodconsumed.y <- ifelse(nutrients_piperine_ttest_total$foodconsumed.y < 0, 0, nutrients_piperine_ttest_total$foodconsumed.y)

result_nutrients_piperine <- t.test(nutrients_piperine_ttest_total$foodconsumed.x, nutrients_piperine_ttest_total$foodconsumed.y, mu = 0, paired = TRUE)
result_nutrients_piperine

nutrients_piperine_ttest_total$difference <- (nutrients_piperine_ttest_total$foodconsumed.x - nutrients_piperine_ttest_total$foodconsumed.y)

# Create a horizontal line plot with colors for positive and negative values
plot_nutrients_piperine <- ggplot(nutrients_piperine_ttest_total, aes(y = difference)) +
  theme_test(base_size = 10) + 
  ggtitle("(C) Trial 1.3") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "bottom") + 
  geom_violin(aes(x = 0, y = difference), fill = "light gray", color = "white", alpha = 0.4) + 
  geom_point(aes(x = 0, color = factor(difference > 0)), size = 5, alpha = 0.7) +
  annotate("segment", x = -Inf, xend = Inf, y = 0, yend = 0, color = "black", linetype = "dashed") + 
  xlab (" ") + 
  ylab ("High nutrients, 2% piperine − Low nutrients, 0.1% piperine") + 
  scale_x_continuous(breaks = NULL) +  # Remove x-axis tick marks
  guides(color = guide_legend(title = NULL)) +  # Remove legend title
  scale_color_manual(values = c("TRUE" = "#440154FF", "FALSE" = "darkorange"),
                     labels = c("Low nutrients,\n0.1% piperine\npreferred", "High nutrients,\n2% piperine\npreferred"),
                     name = "Preference") + 
  scale_y_continuous(breaks = c(-5, 0, 5),
                     limits = c(-5, 5))  # Set the limits of the Y-axis
plot_nutrients_piperine

### Add all graphs together

figure1 <- ggarrange(plot_nutrients,
                     plot_piperine,
                     plot_nutrients_piperine,
                     ncol = 3, nrow = 1,
                     align = "hv")
figure1
ggsave(file="figure1.jpg", 
       plot=figure1,
       width=8,height=5,units="in",dpi=600)

### calculate the amount of nutrients (protein and carbohydrates).
head(nutrients_piperine_ttest_total)

nutrients_piperine_ttest_total$nutrients_dietX <- ((nutrients_piperine_ttest_total$foodconsumed.x)*179/5000)*1000
nutrients_piperine_ttest_total$nutrients_dietY <- ((nutrients_piperine_ttest_total$foodconsumed.y)*65/5000)*1000
nutrients_piperine_ttest_total$nutrient_total <- nutrients_piperine_ttest_total$nutrients_dietX + nutrients_piperine_ttest_total$nutrients_dietY
t1.3_nutrients <- c(nutrients_piperine_ttest_total$nutrient_total)

nutrients_piperine_ttest_total$piperine_dietX <- ((nutrients_piperine_ttest_total$foodconsumed.x)*14.5/5000)*1000
nutrients_piperine_ttest_total$piperine_dietY <- ((nutrients_piperine_ttest_total$foodconsumed.y)*0.7/5000)*1000
nutrients_piperine_ttest_total$piperine_total <- nutrients_piperine_ttest_total$piperine_dietX + nutrients_piperine_ttest_total$piperine_dietY
t1.3_piperine <- c(nutrients_piperine_ttest_total$piperine_total)

###############################
###############################

# create new df with nutrients and piperine consumed per treatment 

head(nutrients_ttest_total)
nutrients_geometric <- nutrients_ttest_total[c("bat", "nutrient_total")]
nutrients_geometric <- nutrients_geometric %>%
  rename(nutrient_total_t1.1 = nutrient_total)

head(piperine_ttest_total)
piperine_geometric <- piperine_ttest_total[c("bat", "nutrient_total", "piperine_total")]
piperine_geometric <- piperine_geometric %>%
  rename(nutrient_total_t1.2 = nutrient_total,
         piperine_total_t1.2 = piperine_total)

head(nutrients_piperine_ttest_total)
nutrients_piperine_geometric <- nutrients_piperine_ttest_total[c("bat", "nutrient_total", "piperine_total")]
nutrients_piperine_geometric <- nutrients_piperine_geometric %>%
  rename(nutrient_total_t1.3 = nutrient_total,
         piperine_total_t1.3 = piperine_total)

geometric <- nutrients_geometric %>%
  left_join(piperine_geometric, by = "bat") %>%
  left_join(nutrients_piperine_geometric, by = "bat") 

head(geometric)

###

# Create a data frame for plotting
plot_data <- rbind(
  data.frame(x = geometric$nutrient_total_t1.1, y = 0, Group = "Trial 1.1"),
  data.frame(x = geometric$nutrient_total_t1.2, y = geometric$piperine_total_t1.2, Group = "Trial 1.2"),
  data.frame(x = geometric$nutrient_total_t1.3, y = geometric$piperine_total_t1.3, Group = "Trial 1.3")
)

plot_data <- rbind(
  data.frame(x = geometric$nutrient_total_t1.1, y = 0, Group = "Trial 1.1", bat = geometric$bat),
  data.frame(x = geometric$nutrient_total_t1.2, y = geometric$piperine_total_t1.2, Group = "Trial 1.2", bat = geometric$bat),
  data.frame(x = geometric$nutrient_total_t1.3, y = geometric$piperine_total_t1.3, Group = "Trial 1.3", bat = geometric$bat)
)

# Calculate mean and standard deviation per group
summary_data <- plot_data %>%
  group_by(Group) %>%
  summarize(Mean_x = mean(x, na.rm = TRUE), 
            Mean_y = mean(y, na.rm = TRUE),
            SD_x = sd(x, na.rm = TRUE),
            SD_y = sd(y, na.rm = TRUE))

# Create the scatterplot
geometric_plot <- ggplot() +
  theme_classic(base_size = 15) + 
  geom_point(data = plot_data, aes(x = x, y = y, color = Group), size = 2.5, alpha = 0.5) +
  geom_point(data = summary_data, aes(x = Mean_x, y = Mean_y, color = Group), size = 4) +
  geom_errorbar(data = summary_data, aes(x = Mean_x, ymin = Mean_y - SD_y, ymax = Mean_y + SD_y, color = Group), 
                width = 9) +
  geom_errorbarh(data = summary_data, aes(y = Mean_y, xmin = Mean_x - SD_x, xmax = Mean_x + SD_x, color = Group), 
                 height = 0.3) +
  labs(x = "Nutrients consumed per trial (mg)",
       y = "Piperine consumed per trial (mg)") +
  scale_color_manual(values = c("Trial 1.1" = "darkolivegreen3", "Trial 1.2" = "#FDE725FF", "Trial 1.3" = "#440154FF")) +
  guides(color = guide_legend(title = "Trials")) 

geometric_plot

data_geometric_nutrients <- plot_data %>%
  filter(Group %in% c("Trial 1.1"))

summary_data <- data_geometric_nutrients %>%
  group_by(Group) %>%
  summarize(Mean_x = mean(x, na.rm = TRUE), 
            Mean_y = mean(y, na.rm = TRUE),
            SD_x = sd(x, na.rm = TRUE),
            SD_y = sd(y, na.rm = TRUE))

geometric_plot_nutrients <- ggplot() +
  theme_test(base_size = 10) + 
  ggtitle("(A) Trial 1.1")+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_point(data = data_geometric_nutrients, aes(x = x, y = y), color = "darkolivegreen3", size = 5, alpha = 0.7) +
  geom_errorbar(data = summary_data, aes(x = Mean_x, ymin = Mean_y - SD_y, ymax = Mean_y + SD_y), 
                color = "darkolivegreen3", width = 9) +
  geom_errorbarh(data = summary_data, aes(y = Mean_y, xmin = Mean_x - SD_x, xmax = Mean_x + SD_x), 
                 color = "darkolivegreen3", height = 0.3) +
  labs(x = "Nutrients consumed (mg)",
       y = "Piperine consumed (mg)") +
  geom_point(data = summary_data, aes(x = Mean_x, y = Mean_y), color = "darkolivegreen4", size = 6, shape = 17) +
  scale_y_continuous(breaks = c(0, 5, 10, 15),
                     limits = c(-1, 17)) + 
  scale_x_continuous(breaks = c(0, 150, 300),
                     limits = c(0, 300)) + 
  geom_text(data = NULL, aes(x = 179, y = 0, label = 'High nutrient\noption'), size = 3, vjust = -0.5) +
  geom_text(data = NULL, aes(x = 65, y = 0 , label = 'Low nutrient\noption'), size = 3, vjust = -0.5) + 
  geom_point(data = NULL, aes(x = 179, y = 0), color = 'black', size = 4) +
  geom_point(data = NULL, aes(x = 65, y = 0), color = 'black', size = 4)

geometric_plot_nutrients

###
data_geometric_piperine <- plot_data %>%
  filter(Group %in% c("Trial 1.2"))

summary_data <- data_geometric_piperine %>%
  group_by(Group) %>%
  summarize(Mean_x = mean(x, na.rm = TRUE), 
            Mean_y = mean(y, na.rm = TRUE),
            SD_x = sd(x, na.rm = TRUE),
            SD_y = sd(y, na.rm = TRUE))

geometric_plot_piperine <- ggplot() +
  theme_test(base_size = 10) + 
  ggtitle("(B) Trial 1.2")+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_point(data = data_geometric_piperine, aes(x = x, y = y), color = "goldenrod1", size = 5, alpha = 0.7) +
  geom_errorbar(data = summary_data, aes(x = Mean_x, ymin = Mean_y - SD_y, ymax = Mean_y + SD_y), 
                color = "goldenrod1", width = 9) +
  geom_errorbarh(data = summary_data, aes(y = Mean_y, xmin = Mean_x - SD_x, xmax = Mean_x + SD_x), 
                 color = "goldenrod1", height = 0.3) +
  labs(x = "Nutrients consumed (mg)",
       y = " ") +
  geom_point(data = summary_data, aes(x = Mean_x, y = Mean_y), color = "goldenrod2", size = 6, shape = 17) + 
  scale_y_continuous(breaks = c(0, 5, 10, 15),
                     limits = c(-1, 17)) + 
  scale_x_continuous(breaks = c(0, 150, 300),
                     limits = c(-1, 300)) + 
  geom_text(data = NULL, aes(x = 245, y = 14.5, label = '2% piperine\n option'), size = 3) +
  geom_text(data = NULL, aes(x = 245, y = 0.7, label = '0.1% piperine\n option'), size = 3) + 
  geom_segment(data = NULL, aes(x = 179, xend = 179, y = 0, yend = 14.5), linetype = 'dashed', color = 'grey') +
  geom_segment(data = NULL, aes(x = 0, xend = 179, y = 14.5, yend = 14.5), linetype = 'dashed', color = 'grey') +
  geom_segment(data = NULL, aes(x = 179, xend = 179, y = 0, yend = 0.7), linetype = 'dashed', color = 'grey') +
  geom_segment(data = NULL, aes(x = 0 , xend = 179, y = 0.7, yend = 0.7), linetype = 'dashed', color = 'grey') + 
  geom_point(data = NULL, aes(x = 179, y = 14.5), color = 'black', size = 4) +
  geom_point(data = NULL, aes(x = 179, y = 0.7), color = 'black', size = 4)

geometric_plot_piperine

###

data_geometric_nutrientspiperine <- plot_data %>%
  filter(Group %in% c("Trial 1.3"))

summary_data <- data_geometric_nutrientspiperine %>%
  group_by(Group) %>%
  summarize(Mean_x = mean(x, na.rm = TRUE), 
            Mean_y = mean(y, na.rm = TRUE),
            SD_x = sd(x, na.rm = TRUE),
            SD_y = sd(y, na.rm = TRUE))

geometric_plot_nutrientspiperine <- ggplot() +
  theme_test(base_size = 10) + 
  ggtitle("(C) Trial 1.3")+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_point(data = data_geometric_nutrientspiperine, aes(x = x, y = y), color = "pink3", size = 5, alpha = 0.7) +
  geom_errorbar(data = summary_data, aes(x = Mean_x, ymin = Mean_y - SD_y, ymax = Mean_y + SD_y), 
                color = "pink3", width = 9) +
  geom_errorbarh(data = summary_data, aes(y = Mean_y, xmin = Mean_x - SD_x, xmax = Mean_x + SD_x), 
                 color = "pink3", height = 0.3) +
  labs(x = "Nutrients consumed (mg)",
       y = " ") +
  geom_point(data = summary_data, aes(x = Mean_x, y = Mean_y), color = "pink4", size = 6, shape = 17) +
  scale_y_continuous(breaks = c(0, 5, 10, 15),
                     limits = c(-1, 17)) + 
  scale_x_continuous(breaks = c(0, 150, 300),
                     limits = c(0, 300)) + 
  geom_text(data = NULL, aes(x = 179, y = 14.5, label = 'High nutrient,\n2% piperine option'), size = 3, vjust = -0.5) +
  geom_text(data = NULL, aes(x = 80, y = 0.7, label = 'Low nutrient,\n0.1% piperine option'), size = 3, vjust = -0.5) + 
  geom_segment(data = NULL, aes(x = 179, xend = 179, y = 0, yend = 14.5), linetype = 'dashed', color = 'grey') +
  geom_segment(data = NULL, aes(x = 0, xend = 179, y = 14.5, yend = 14.5), linetype = 'dashed', color = 'grey') +
  geom_segment(data = NULL, aes(x = 65, xend = 65, y = -0, yend = 0.7), linetype = 'dashed', color = 'grey') +
  geom_segment(data = NULL, aes(x = 0, xend = 65, y = 0.7, yend = 0.7), linetype = 'dashed', color = 'grey') + 
  geom_point(data = NULL, aes(x = 179, y = 14.5), color = 'black', size = 4) +
  geom_point(data = NULL, aes(x = 65, y = 0.7), color = 'black', size = 4)

geometric_plot_nutrientspiperine

figure2 <- ggarrange(geometric_plot_nutrients,
                             geometric_plot_piperine,
                             geometric_plot_nutrientspiperine,
                             ncol = 3, nrow = 1,
                             align = "hv")
figure2 
ggsave(file="figure2.jpg", 
       plot=figure2 ,
       width=8,height=4,units="in",dpi=600)

#LM
head(plot_data)
plot_data_cleaned <- plot_data[complete.cases(plot_data$x, plot_data$y), ]
hist(plot_data$y)
hist(plot_data$x)

geometric_nutrients_gl <- lm(x ~ Group, data = plot_data_cleaned)
summary(geometric_nutrients_gl)
plot(allEffects(geometric_nutrients_gl))
shapiro.test(resid(geometric_nutrients_gl))
emmeans(geometric_nutrients_gl, pairwise ~ Group)
r2(geometric_nutrients_gl)
parameters(geometric_nutrients_gl)

glmm_piperine <- plot_data_cleaned %>%
  filter(!Group %in% c("Trial 1.1"))

geometric_piperine_gl <- lm(y ~ Group, data = glmm_piperine)
summary(geometric_piperine_gl)
plot(allEffects(geometric_piperine_gl))
shapiro.test(resid(geometric_piperine_gl))
parameters(geometric_piperine_gl)
r2(geometric_piperine_gl)

###################################
###################################
###################################
###################################
###################################
###################################
###################################

