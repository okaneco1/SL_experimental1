# Statistical Analyses for Experiment 1 -- Temperature/Fasting Periods

options(scipen = 999)

# libraries
library(COUNT)
library(tidyverse)
library(readxl)
library(pscl)
library(lme4)
library(RColorBrewer)
library(multcomp)
library(ggpubr)


#----- data setup

# load data
full_data <- read_csv("Experiment_1_full_data_table.csv")

# data adjustments - combining S. namaycush and Salmonidae reads
#full_data <- full_data %>%
  #mutate(trout = Salvelinus_namaycush + Salmonidae_unclassified)
colnames(full_data)[20] <- "weight_gain"
# removing "_mean" from column names
colnames(full_data) <- sub("_mean", "", colnames(full_data))
# change column names
colnames(full_data)[colnames(full_data) == "Fasting period(days)"] <- "fasting_period"
colnames(full_data)[colnames(full_data) == "Weight loss(g)"] <- "weight_loss"
colnames(full_data)[colnames(full_data) == "Tube ID"] <- "sample"
# columnb for temp as a factor
full_data$temp_f <- as.factor(full_data$temp)

# using only lake trout and sea lamprey reads
full_data <- full_data[,c(1:2,6:7,14:ncol(full_data))]
# column for total reads in sample
full_data$total_reads <- rowSums(full_data[, 2:4], na.rm = TRUE)



#---- combining host fish reads

# all trout
full_data <- full_data %>%
  mutate(all_trout = Salvelinus_namaycush + Salmonidae_unclassified)



#----- basic data

# initial weight data
mean(full_data$`Initial Weight(g)`) # 6.475714
sd(full_data$`Initial Weight(g)`) # 0.9888862



# -----------------------  STATISTICAL ANALYSES -----------------------------

# aov(lm(...))

#----------- ANOVA and pairwise comparisons -----------

aov1 <- aov(all_trout ~ as.factor(temp), data = full_data) 
aov2 <- aov(all_trout ~ temp_f * as.factor(fasting_period), data = full_data) 
aov3 <- aov(all_trout ~ as.factor(temp) * as.factor(fasting_period) + weight_gain, data = full_data)
#aov2_log <- aov(log(all_trout) ~ as.factor(temp) * as.factor(fasting_period), data = full_data) 

summary(aov2)

anova(aov2, aov3) # weight gain might add significant improvement to model

#---- testing assumptions
library(car)
#verifying normality
plot(aov2, which = 2)
qqPlot(aov2$residuals, id = FALSE)

hist(resid(aov2), breaks = 10)

# homogeneity of variances
plot(aov2, which = 3) # not great...


#---- post-hoc Tukey's HSD 
library(agricolae)
# full model
TukeyHSD(aov2)
# subsets
TukeyHSD(aov2, which = "as.factor(fasting_period)")
TukeyHSD(aov2, which = "temp_f:as.factor(fasting_period)")

# only sig dif comparisons from interactions
tukey_results <- TukeyHSD(aov2)
interaction_results <- tukey_results$`temp_f:as.factor(fasting_period)`
# convert to dataframe
int_df <- as.data.frame(interaction_results)
int_df_sig <- filter(int_df, `p adj` < 0.05)
#ordered by diff
int_df_sig_ordered <- int_df_sig[order(int_df_sig$diff), ]
#--- write out
#int_df_sig_ordered$comparisons <- row.names(int_df_sig_ordered)
#write_csv(int_df_sig_ordered, "sig_dif_interactions.csv")

# denotes letter for groups that are not sig dif from each other
HSD.test(aov2, trt = c("as.factor(temp)", "as.factor(fasting_period)"), console = TRUE)
# plot
par(mar = c(4.5, 7, 3, 2))
#plot(TukeyHSD(aov2, which = "as.factor(temp):as.factor(fasting_period)"), las = 2)




#--------------- Linear Regressions ---------

# weight change models

# weight gain and sequence read count
lm_wg <- lm(all_trout/total_reads ~ weight_gain, data = full_data)
summary(lm_wg)

# weight loss and fasting period
lm_wl_fp <- lm(weight_loss ~ as.factor(fasting_period), data = full_data)
summary(lm_wl_fp) # 0.0278
anova(lm_wg_temp)
TukeyHSD(aov(lm(weight_loss ~ as.factor(fasting_period), data = full_data)))

# weight gain and temperature
lm_wg_temp <- lm(weight_gain ~ as.factor(temp), data = full_data)
summary(lm_wg_temp)
anova(lm_wg_temp) # 0.0000000001517
TukeyHSD(aov(lm(weight_gain ~ as.factor(temp), data = full_data)))



#----------------- VISUALIZATIONS ------------------------

# Full Comparison Between Temp and Fasting Period

# viewing data
sorted_data <- full_data %>%
  dplyr::select(sample, all_trout, Petromyzontidae_unclassified, temp, fasting_period) %>%
  dplyr::arrange(temp, fasting_period)


# first split by temp
custom_labels <- c('5' = 'Temperature = 5°C', '10' = 'Temperature = 10°C', '15' = 'Temperature = 15°C')

ggplot(full_data, aes(x = fasting_period, y = all_trout/total_reads, color = factor(temp))) +
  geom_boxplot(aes(group = interaction(fasting_period, factor(temp))),
               alpha = 0.5, position = position_dodge(width = 3), 
               outlier.shape = NA,  coef = Inf) +
  geom_smooth(method = "lm", aes(group = factor(temp)), formula = y ~ x, se = FALSE) +
  labs(color = "Temperature (°C)", y = "Total Lake Trout Reads", x = "Fasting Period (Days)") +
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_brewer(palette = "Dark2") +
  scale_x_continuous(breaks = c(0, 5, 10, 20, 30)) +
  facet_wrap(vars(temp), labeller = labeller(temp = custom_labels)) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold", size = 12))

# then split by fasting period
custom_labels2 <- c('0' = '0 Days', '5' = '5 Days', '10' = '10 Days', '20' = '20 Days', '30' = '30 Days')

ggplot(full_data, aes(x = temp, y = all_trout, fill = factor(fasting_period))) +
  geom_boxplot(aes(group = interaction(temp, factor(fasting_period))),
               alpha = 0.5, position = position_dodge(width = 3), 
               outlier.shape = NA, coef = Inf) +
  #geom_smooth(method = "lm", aes(group = factor(fasting_period)), formula = y ~ x, se = FALSE) +
  labs(fill = "Fasting Period (days)", y = "Sequence Read Count (Lake Trout)", x = "Temperature (°C)") +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_brewer(palette = "Dark2") +
  #scale_x_continuous(breaks = c(0, 5, 10, 20, 30)) +
  facet_wrap(~fasting_period, nrow = 1,scales = "free_x", labeller = labeller(fasting_period = custom_labels2)) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold", size = 12))

# another way to plot
ggline(full_data, x = "fasting_period", y = "all_trout", color = "temp_f", add = c("mean_se")) +
  labs(y = "Lake Trout Sequence Read Count", 
       x = "Fasting Period (Days)",
       color = "Temperature (°C)") +
  scale_color_brewer(palette = "Set2")



 
# weight gain comparison
cor_coeff_wg <- cor(full_data$weight_gain, full_data$Salvelinus_namaycush, use = "complete.obs")

ggplot(full_data, aes(x = weight_gain, y = Salvelinus_namaycush)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "deepskyblue") +
  labs(x = "Weight Gain (g)", y = "Lake Trout Sequence Read Count") +
  theme_minimal() +
  annotate("text", x = Inf, y = Inf, label = paste("R =", round(cor_coeff_wg, 3)), 
           hjust = 2, vjust = 15, size = 5, color = "deepskyblue")
# does not seem to be relationship between weight gain and sequence read count

# weight loss and fasting period
cor_coeff_wl <- cor(full_data$`Weight loss(g)`, full_data$Salvelinus_namaycush, use = "complete.obs")

ggplot(full_data, aes(x = fasting_period, y = `Weight loss(g)`)) +
  geom_boxplot(aes(group = fasting_period), outlier.shape = NA) +
  #geom_smooth(method = "lm", se = FALSE, color = "deepskyblue") +
  geom_jitter(width = 0.2, alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgreen", alpha = 0.5) +
  theme_minimal() 


# weight gain and temperature
ggplot(full_data, aes(x = temp, y = weight_gain)) +
  geom_boxplot(aes(group = temp), outlier.shape = NA) +
  #geom_smooth(method = "lm", se = FALSE, color = "deepskyblue") +
  geom_jitter(width = 0.2, alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgreen", alpha = 0.5) +
  labs(x = "Temperature (°C)", y = "Weight Gain (g)") +
  theme_minimal()


# more in-depth plot of Tukey HSD interaction results
tukey_result <- TukeyHSD(aov2, "as.factor(temp):as.factor(fasting_period)")
tukey_df <- as.data.frame(tukey_result$`as.factor(temp):as.factor(fasting_period)`)
tukey_df$comparison <- rownames(tukey_df)
# identify non-zero crossing intervals
tukey_df$highlight <- ifelse(tukey_df$lwr > 0 | tukey_df$upr < 0, "Significant", "Non-Significant")
# plot
ggplot(tukey_df, aes(x = comparison, y = diff)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr, ymax = upr, color = highlight), width = 0.2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip() +
  labs(color = "") +
  geom_hline(yintercept = 0, color = "black", alpha = 0.5) +
  scale_color_manual(values = c("Significant" = "red", "Non-Significant" = "black")) 
