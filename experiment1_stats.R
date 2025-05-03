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
library(MASS)
library(pscl)


#----- data setup

# load data
full_data <- read_csv("Experiment_1_full_data_table.csv")

# data adjustments - combining S. namaycush and Salmonidae reads
#full_data <- full_data %>%
  #mutate(trout = Salvelinus_namaycush + Salmonidae_unclassified)
colnames(full_data)[21] <- "weight_gain"
# removing "_mean" from column names
colnames(full_data) <- sub("_mean", "", colnames(full_data))
# change column names
colnames(full_data)[colnames(full_data) == "Fasting period(days)"] <- "fasting_period"
colnames(full_data)[colnames(full_data) == "Weight loss(g)"] <- "weight_loss"
colnames(full_data)[colnames(full_data) == "Tube ID"] <- "sample"
colnames(full_data)[colnames(full_data) == "Days Attached"] <- "days_attached"
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

# organize into table for visualization
lt_sequences_table_simplified <- full_data %>%
  select(sample, all_trout, temp, fasting_period) %>%
  arrange(temp, fasting_period)

# set labels for fasting days
fasting_labels <- lt_sequences_table_simplified %>%
  distinct(sample, fasting_period) %>%
  deframe() 

# refactor fasting days for proper order
lt_sequences_table_simplified <- lt_sequences_table_simplified %>%
  mutate(fasting_period = as.numeric(fasting_period)) %>%
  mutate(sample = factor(sample, levels = unique(sample[order(fasting_period)])))

# visualization
ggplot(lt_sequences_table_simplified, aes(x = sample, y = all_trout, fill = as.factor(temp)))+
  geom_col() +
  scale_x_discrete(labels = fasting_labels) +
  labs(x = "Fasting Days",
       y = "Lake Trout Sequence Reads") +
  theme_minimal(base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,),
        panel.grid.major = element_line(color = alpha("gray", 0.1)),
        panel.grid.minor = element_line(color = alpha("gray", 0.1))
        )



#----- basic data

# initial weight data
mean(full_data$`Initial Weight(g)`) # 6.475714
sd(full_data$`Initial Weight(g)`) # 0.9888862



# -----------------------  STATISTICAL ANALYSES -----------------------------
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
anova(lm_wl_fp)
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
ggline(full_data, x = "fasting_period", y = "all_trout", 
       color = "temp_f", add = c("mean_se"), position = position_dodge(width = 0.1)) +
  labs(y = "Lake Trout Sequence Read Count", 
       x = "Fasting Period (Days)",
       #title = "Lake Trout Sequence Read Count vs Fasting Days",
       color = "Temperature (°C)") +
  theme(legend.position = c(0.7, 0.8),
        plot.title = element_text(face = "bold"),  
        axis.title.x = element_text(face = "bold"),  
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_color_brewer(palette = "Set2")

# checking out outliers in 15°C / 30 days fasting
ggplot(full_data, aes(x = fasting_period, y = all_trout, color = temp_f)) +
  geom_jitter(width = 0.4, height = 0.1)

filter(full_data, temp == 15, )
full_data[full_data$fasting_period == 30 & full_data$temp_f == 15 & full_data$all_trout > 5000,]

 
# weight gain comparison
cor.test(full_data$weight_gain, full_data$Salvelinus_namaycush, use = "complete.obs")
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
cor_coeff_wl <- cor(full_data$weight_loss, full_data$Salvelinus_namaycush, use = "complete.obs")

ggplot(full_data, aes(x = fasting_period, y = weight_loss)) +
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
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16), 
    axis.title.x = element_text(face = "bold"),          
    axis.title.y = element_text(face = "bold")           
  )


# weight and days attached
cor.test(full_data$days_attached, full_data$weight_gain)
correlation <- cor(full_data$days_attached, full_data$weight_gain)

ggplot(full_data, aes(x = days_attached, y = weight_gain)) +
  geom_point(aes(group = temp)) +
  geom_smooth(method = "lm", se = FALSE, color = "deepskyblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgreen", alpha = 0.5) +
  labs(x = "Days Attached", y = "Weight Gain (g)") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16), 
    axis.title.x = element_text(face = "bold"),          
    axis.title.y = element_text(face = "bold")           
  ) #+
  #annotate("text", x = max(full_data$days_attached) - 5, y = max(full_data$weight_gain) - 1, 
           #label = paste("Correlation: ", round(correlation, 2)), 
           #color = "black", size = 4, hjust = 1)



temp_weight_anova <- aov(weight_gain ~ as.factor(temp), data = full_data)
summary(temp_weight_anova)

tukey_results <- TukeyHSD(temp_weight_anova)
plot(tukey_results)


# ------- MODEL SELECTION

# viewing data
hist(log(full_data$all_trout))
hist(full_data$all_trout)

# plot each variable against read count 
ggplot(full_data, aes(x = weight_gain, y = all_trout)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()
ggplot(full_data, aes(x = fasting_period, y = all_trout)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, span = 1) + # span increased to desensitize fit
  labs(x = "Fasting Period (Days)", y = "Lake Trout Read Count",
       title = "Lake Trout Sequence Reads per Fasting Period") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14), 
    axis.title.x = element_text(face = "bold", size = 13),          
    axis.title.y = element_text(face = "bold", size = 13)           
  )
ggplot(full_data, aes(x = temp, y = all_trout)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()


# setting up  models
lm1 <- lm(all_trout ~ weight_gain + fasting_period * temp + days_attached, data = full_data)
poisson_model <- glm(all_trout ~ weight_gain + fasting_period * temp + days_attached, family = poisson, data = full_data)
negbin_lm <- glm.nb(round(all_trout) ~ weight_gain + temp * fasting_period + days_attached, data = full_data)
zinb_lm <- zeroinfl(round(all_trout) ~ weight_gain + fasting_period * as.factor(temp) + days_attached | 1, 
                    data = full_data, 
                    dist = "negbin")

# summaries
summary(lm1)
summary(poisson_model)
summary(negbin_lm)
summary(zinb_lm)

# AIC/model comparisons
AIC(lm1)
AIC(poisson_model)
AIC(negbin_lm)
AIC(zinb_lm)

vuong(zinb_lm, negbin_lm)

# as there are not many zeros in the averaged read count data (as compared to 
# replicate-specific data), the negative binomial model is preferred over the
# zero-inflated negative binomial

(1 - exp(-0.021))*100 # 2.078104

# visual table for model comparisons
model_names <- c("Linear Model", "Poisson GLM", "Negative Binomial GLM", "Zero-Inflated NB")
aic_values <- c(
  AIC(lm1),
  AIC(poisson_model),
  AIC(negbin_lm),
  AIC(zinb_lm)
)

# create data frame
aic_table <- data.frame(
  Model = model_names,
  AIC = aic_values
)

# calculate ΔAIC and AIC weights
aic_table$Delta_AIC <- aic_table$AIC - min(aic_table$AIC)
aic_table$AIC_weight <- exp(-0.5 * aic_table$Delta_AIC)
aic_table$AIC_weight <- aic_table$AIC_weight / sum(aic_table$AIC_weight)

# sort 
aic_table <- aic_table[order(aic_table$AIC), ]

# round
aic_table$AIC <- round(aic_table$AIC, 2)
aic_table$Delta_AIC <- round(aic_table$Delta_AIC, 2)
aic_table$AIC_weight <- round(aic_table$AIC_weight, 3)

# table and output
aic_table
write.csv(aic_table, "aic_model_comparison.csv", row.names = FALSE)


full_data %>%
  group_by(temp, fasting_period) %>%
  summarize(
    mean_weight_gain = mean(weight_gain, na.rm = TRUE),
    sd = sd(weight_gain, na.rm = TRUE),
    sample_size = n(),  # This adds the sample size for each group
    .groups = 'drop'
  )
mean(full_data$weight_gain)



# can check different variable combinations for the neg bin model
negbin_lm2 <- glm.nb(round(all_trout) ~ weight_gain, data = full_data)
negbin_lm3 <- glm.nb(round(all_trout) ~ temp, data = full_data)
negbin_lm4 <- glm.nb(round(all_trout) ~ fasting_period, data = full_data)
negbin_lm5 <- glm.nb(round(all_trout) ~ days_attached, data = full_data)
negbin_lm6 <- glm.nb(round(all_trout) ~ temp * fasting_period, data = full_data)
negbin_lm7 <- glm.nb(round(all_trout) ~ weight_gain + temp * fasting_period, data = full_data)
negbin_lm8 <- glm.nb(round(all_trout) ~ temp * fasting_period + days_attached, data = full_data)
negbin_lm9 <- glm.nb(round(all_trout) ~ weight_gain + temp * fasting_period + days_attached, data = full_data)
negbin_lm_null <- glm.nb(round(all_trout) ~ 1, data = full_data)

aic_df <- AIC(negbin_lm, negbin_lm2, negbin_lm3, negbin_lm4, negbin_lm5, 
    negbin_lm6, negbin_lm7, negbin_lm8, negbin_lm9, negbin_lm_null)

model_descriptions <- c("weight_gain + temp + fasting_period + days_attached",
                        "weight_gain",
                        "temp",
                        "fasting_period",
                        "days_attached",
                        "temp * fasting_period",
                        "weight_gain + temp * fasting_period",
                        "temp * fasting_period + days_attached",
                        "weight_gain + temp * fasting_period + days_attached",
                        "null_model")

aic_df$model <- model_descriptions

arrange(aic_df, AIC)

# add delta AIC
aic_df <- aic_df %>%
  arrange(AIC) %>%
  mutate(Delta_AIC = round(AIC - min(AIC), 2))
aic_df_reorder <- aic_df[c(3,1,2,4)]

# write out
write.csv(aic_df_reorder, file = "experiment1_negbin_model_selection.csv")



# coefficients
best_model <- glm.nb(round(all_trout) ~ temp * fasting_period, data = full_data)
summary(negbin_lm6)

1-exp(-0.118774) #temp
1-exp(0.007693) #interaction


# deviance
null_deviance <- best_model$null.deviance
residual_deviance <- best_model$deviance

# calculate explained deviance (pseudo R-squared)
explained_deviance <- 1 - (residual_deviance / null_deviance)
explained_deviance

#mcfaddens R-squared
pR2(best_model)


# using MuMIn and dredge
library(MuMIn)

global_model <- glm.nb(round(all_trout) ~ weight_gain + temp * fasting_period + days_attached, na.action = na.fail, data = full_data)
dredge_model <- dredge(global_model)

dredge_df <- data.frame(dredge_model)
write.csv(dredge_df, file = "negative binomial predictor table exp1.csv", row.names = FALSE)



#------- misc.

# investigating read counts at 30 days
filtered_30_days <- full_data %>%
  filter(fasting_period == 30) %>%
  select(sample, all_trout, fasting_period, temp, )





