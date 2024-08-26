# Categorize Detections Between Replicates

# The objective here is to determine which replicates detected lake trout
# sequences, and to compare whether double, single, or no detections
# correlates with differences in temperature and fasting period

#----- libraries
library(tidyverse)
library(readxl)
library(nnet)
library(broom)
library(patchwork)
library(sjPlot)
library(brglm2)
library(lmtest)
library(geomtextpath)



#----- load replicate data sets
com_mat <- read_excel("submission1A_community_matrix.xls", range = "A1:AW506")
com_mat_r <- read_excel("submission1B_community_matrix.xls")



#----- organize data
colnames(com_mat)[1] <- "sample"
colnames(com_mat_r)[1] <- "sample"

# making all names the same
com_mat$sample <- sub("\\.12S", "", com_mat$sample)
com_mat_r$sample <- sub("_R\\.12S", "", com_mat_r$sample)

# remove Homo sapiens
com_mat <- com_mat[, -4]
com_mat_r <- com_mat_r[, -4]

# add in a column for lake trout (including Salmonidae)
com_mat <- mutate(com_mat, lake_trout = Salvelinus_namaycush + Salmonidae_unclassified)
com_mat_r <- mutate(com_mat_r, lake_trout = Salvelinus_namaycush + Salmonidae_unclassified)

# add in a column for total reads
com_mat$total_reads <- rowSums(com_mat[, c(2:48)], na.rm = TRUE)
com_mat_r$total_reads <- rowSums(com_mat_r[, c(2:48)], na.rm = TRUE)

# simplify to only lake trout column
com_mat <- com_mat[, c(1, 49:50)]
com_mat_r <- com_mat_r[, c(1, 49:50)]

# assign percentage of lake trout reads to total reads
com_mat$lt_rra <- round((com_mat$lake_trout/com_mat$total_reads)*100, 2)
com_mat_r$lt_rra <- round((com_mat_r$lake_trout/com_mat_r$total_reads)*100, 2)

# filter to only experimental samples (keep T23_30 as is only T23 sample)
com_mat <- com_mat %>%
  filter(grepl("^T22", sample) | sample == "T23_30")
com_mat_r <- com_mat_r %>%
  filter(grepl("^T22", sample) | sample == "T23_30")

# create detection column with 20 read detection limit
threshold <- 10 # makes threshold adjustable
rra <- 1

com_mat <- com_mat %>%
  mutate(detection = ifelse(lake_trout <= threshold | lt_rra < 1, 0, 1))
com_mat_r <- com_mat_r %>%
  mutate(detection = ifelse(lake_trout <= threshold | lt_rra < 1, 0, 1))

# function to unify sample names
unify_names <- function(name) {
  # add 0 before single digits; replace underscore with hyphen
  name <- sub("[_-](\\d)$", "-0\\1", name)
  # replace _ with -
  name <- sub("_", "-", name)
  # insert 20 after T
  name <- sub("T", "T20", name)
  return(name)
}

com_mat$sample <- sapply(com_mat$sample, unify_names)
com_mat_r$sample <- sapply(com_mat_r$sample, unify_names)



#------ create data set that merges replicate detection values

# combine columns, indicating replicates
merged_detections <- left_join(com_mat, com_mat_r, by = "sample", suffix = c("_rep1", "_rep2"))
merged_detections[is.na(merged_detections)] <- 0 # change all NAs to 0

# add combined detection column
merged_detections$combined_detection <- NA


# add value for combined detections
for (i in 1:nrow(merged_detections)) {
  # assign detections for replicates of same sample
  rep1_det <- merged_detections$detection_rep1[i]
  rep2_det <- merged_detections$detection_rep2[i]
  
  if (rep1_det + rep2_det == 2) {
    merged_detections$combined_detection[i] <- 2
  } else if (rep1_det + rep2_det == 1) {
    merged_detections$combined_detection[i] <- 1
  } else {
    merged_detections$combined_detection[i] <- 0
  }
}



#------ combine detection dataset with variable data
full_data <- read_csv("Experiment_1_full_data_table.csv")
full_data <- full_data[,c(1,14:25)]
colnames(full_data)[1] <- "sample"

# not all samples from full data are in experimental set
full_detection_data <- inner_join(merged_detections, full_data, by = "sample")



full_detection_data <- full_detection_data %>%
  rename(fasting_period = `Fasting period(days)`,
         weight_loss = `Weight loss(g)`,
         weight_gain = `Weight gain(g)`,
         days_attached = `Days Attached`,
         initial_weight = `Initial Weight(g)`)





#----- statistical analyses

# set categorical variables as factors 
full_detection_data$temp <- factor(full_detection_data$temp,
                                   levels = c(5, 10, 15))
full_detection_data$fasting_period <- factor(full_detection_data$fasting_period,
                                             levels = c(0, 5, 10, 20, 30))
full_detection_data$combined_detection <- factor(full_detection_data$combined_detection,
                                                 levels = c(0, 1, 2))

# trying with binomial distributions
# first, for detection in both replicates

# create additional data column for both detections (and label)
full_detection_data$both_det <- ifelse(full_detection_data$combined_detection == 2, 1, 0)
full_detection_data$both_det <- factor(full_detection_data$both_det,
                                       levels = c(0, 1),
                                       labels = c("no detection", "detection"))


# model selection
# null
m_null <- glm(both_det ~ 1, family = "binomial", data = full_detection_data)

# initial model
m1 <- glm(both_det ~ temp + fasting_period + weight_gain/initial_weight + days_attached, family = "binomial", data = full_detection_data)
m2 <- glm(both_det ~ temp * fasting_period + weight_gain/initial_weight + days_attached, family = "binomial", data = full_detection_data)
m3 <- glm(both_det ~ temp * fasting_period + weight_gain * days_attached, family = "binomial", data = full_detection_data)

anova(m1, m2, test = "LRT") # significant
anova(m2, m3, test = "LRT") 

# compare against simplified models
m_temp <- glm(both_det ~ temp, family = "binomial", data = full_detection_data)
m_fast <- glm(both_det ~ fasting_period, family = "binomial", data = full_detection_data)
m_weight <- glm(both_det ~ weight_gain/initial_weight, family = "binomial", data = full_detection_data)
m_attached <- glm(both_det ~ days_attached, family = "binomial", data = full_detection_data)

anova(m_temp, m2, test = "LRT")
anova(m_fast, m2, test = "LRT")
anova(m_weight, m2, test = "LRT")
anova(m_attached, m2, test = "LRT") # all are significant

# compare against slightly more complex models
m_weight2 <- glm(both_det ~ temp * fasting_period + days_attached, family = "binomial", data = full_detection_data)
m_attached2 <- glm(both_det ~ temp * fasting_period + weight_gain/initial_weight, family = "binomial", data = full_detection_data)

anova(m_weight2, m2, test = "LRT") # weight doesn't seem to improve model
anova(m_attached2, m2, test = "LRT") # significant

# isolate weight compared to temp/fasting period
m_temp_fast <- glm(both_det ~ temp * fasting_period, family = "binomial", data = full_detection_data)
anova(m_temp_fast, m_attached2, test = "LRT") # addition of weight is not significant

# isolate days attached
m_attached3 <- glm(both_det ~ temp * fasting_period + days_attached, family = "binomial", data = full_detection_data)
anova(m_temp_fast, m_attached3, test = "LRT") # addition of days attached is significant

# can set as final model
m_final <- glm(both_det ~ temp * fasting_period + days_attached, family = "binomial", data = full_detection_data)

# is interaction still significant?
m_final_nointer <- glm(both_det ~ temp + fasting_period + days_attached, family = "binomial", data = full_detection_data)
anova(m_final_nointer, m_final, test = "LRT") # yes

# summary
summary(m_final)


# only days_attached is significant here, but this is comparing variable levels of both
# fasting period and temperature, not overall effects of fasting period and temperature

# can use no-interaction model to compare overall effects with isolated effect model
m_final_temp <- glm(both_det ~ fasting_period + days_attached, family = "binomial", data = full_detection_data)
m_final_fast <- glm(both_det ~ temp + days_attached, family = "binomial", data = full_detection_data)

# is temperature overall significant?
anova(m_final_temp, m_final_nointer, test = "LRT") # no
anova(m_final_temp, m_final, test = "LRT") # but, yes if included as an interaction

# is fasting period overall significant?
anova(m_final_fast, m_final_nointer, test = "LRT") # yes
anova(m_final_fast, m_final, test = "LRT") # yes


# so, overall, temperature is significant overall in its interaction with fasting
# period (in detection for both replicates) and fasting period is also significant







#------ can now look for patterns where detection is in *at least one* replicate
# create new columns (and label)
full_detection_data$one_det <- ifelse(full_detection_data$combined_detection == 1 |
                                        full_detection_data$combined_detection == 2, 1, 0)
full_detection_data$one_det <- factor(full_detection_data$one_det,
                                      levels = c(0, 1),
                                      labels = c("no detection", "detection"))

# model selection
m_null <- glm(one_det ~ 1, family = "binomial", data = full_detection_data)

# initial model
m1 <- glm(one_det ~ temp + fasting_period + weight_gain/initial_weight + days_attached, family = "binomial", data = full_detection_data)
m2 <- glm(one_det ~ temp * fasting_period + weight_gain/initial_weight + days_attached, family = "binomial", data = full_detection_data)
m3 <- glm(one_det ~ temp * fasting_period + weight_gain * days_attached, family = "binomial", data = full_detection_data)

anova(m1, m2, test = "LRT") # not significant
anova(m2, m3, test = "LRT") 

summary(m1) # only days attached is significant overall
# can investigate whether temp or fasting period are significant overall

m1_weight <- glm(one_det ~ temp + fasting_period + days_attached, family = "binomial", data = full_detection_data)
m1_temp <- glm(one_det ~ fasting_period + weight_gain/initial_weight + days_attached, family = "binomial", data = full_detection_data)
m1_fasting <- glm(one_det ~ temp + weight_gain/initial_weight + days_attached, family = "binomial", data = full_detection_data)

# are either weight, temp or fasting period significant overall?
anova(m1_weight, m1, test = "LRT") # no
anova(m1_temp, m1, test = "LRT") # no
anova(m1_fasting, m1, test = "LRT") # yes

# only temperature in this case 
# so, try removing weight and fasting from the model and compare via AIC
m1_temp_attach <- glm(one_det ~ temp + days_attached, family = "binomial", data = full_detection_data)

AIC(m1, m1_weight, m1_temp_fast, m1_temp_attach) # temp_attach is best (lowest df and AIC)

# do interactions improve fit?
m1_temp_fast_inter <- glm(one_det ~ temp * fasting_period + days_attached, family = "binomial", data = full_detection_data)
AIC(m1_temp_attach, m1_weight, m1, m1_temp_fast_inter) # no



#-----------------------
# ordinal logistic model selection (with bias reduction)

# null
olm_null <- bracl(combined_detection ~ 1, data = full_detection_data)
# test interaction
olm1 <- bracl(combined_detection ~ temp + fasting_period + days_attached, data = full_detection_data)
olm2 <- bracl(combined_detection ~ temp * fasting_period + days_attached, data = full_detection_data)
# test significance of variables
olm3 <- bracl(combined_detection ~ fasting_period + days_attached, data = full_detection_data)
olm4 <- bracl(combined_detection ~ temp + days_attached, data = full_detection_data)
olm5 <- bracl(combined_detection ~ temp + fasting_period, data = full_detection_data)

# using likelihood ratio test to compare nested models for overall significance
# interaction significance
lrtest(olm1, olm_null) # improvement over null
lrtest(olm2, olm1) # interaction is not significant

# variable significance
lrtest(olm1, olm3) # temp not significant
lrtest(olm1, olm4) # fasting period is significant
lrtest(olm1, olm5) # days attached is significant

summary(olm1)




#------ visualize interactions
temp_labels <- c('5' = 'Temperature = 5°C', '10' = 'Temperature = 10°C', '15' = 'Temperature = 15°C')
fasting_labels <- c('0' = '0 Days Fasted',
                    '5' = '5 Days Fasted',
                    '10' = '10 Days Fasted',
                    '20' = '20 Days Fasted',
                    '30' = '30 Days Fasted')


# faceted temperature
ggplot(full_detection_data, aes(x = factor(fasting_period), fill = factor(combined_detection))) +
  geom_bar(position = "fill", stat = "count") + 
  labs(
    fill = "Detections", 
    y = "Proportion of Lake Trout Detections", 
    x = "Fasting Period (Days)",
    title = "Proportion of Detections between Two Replicates") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(
    breaks = c("0", "5", "10", "20", "30")) +
  facet_wrap(vars(temp), labeller = labeller(temp = temp_labels)) +
  theme(
    strip.text = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 12))

# faceted fasting period
ggplot(full_detection_data, aes(x = factor(temp), fill = factor(combined_detection))) +
  geom_bar(position = "fill", stat = "count") + 
  labs(
    fill = "Detections", 
    y = "Proportion of Lake Trout Detections", 
    x = "Temperature (C°)",
    title = "Proportion of Detections between Two Replicates") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(
    breaks = c("5", "10", "15")) +
  facet_wrap(vars(fasting_period), labeller = labeller(temp = fasting_labels), nrow = 1) +
  theme(
    strip.text = element_text(size = 10),
    plot.title = element_text(face = "bold"), 
    axis.title.x = element_text(face = "bold"), 
    axis.title.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold"))


# combining plots
#combined_plot <- (rra_0 + rra_0.5 + rra_1 + rra_2 + rra_4) + 
#plot_layout(guides = "collect") 




#------ OCCUPANCY MODEL
library(unmarked)

# set up detection matrix (replicates)
detections <- data.frame(rep1 = full_detection_data$detection_rep1,
                         rep2 = full_detection_data$detection_rep2)

# adding covariates, two types:
#--- observation level (weight gain?)
#--- site level (temp and fasting period?)

# set up site covariates (can think of these as don't change between replicates)
site_covs <- as.data.frame(full_detection_data[ ,c(14,18:20)]) # temp and fasting period
#site_covs$temp <- as.numeric(as.character(site_covs$temp))
site_covs$fasting_period <- as.numeric(as.character(site_covs$fasting_period))

# set up observation covariates
#obs_covs <- list(weight_gain = full_detection_data[,18])

# create covariate occupancy frame object
cov_occu <- unmarkedFrameOccu(y = detections, siteCovs = site_covs)
summary(cov_occu)

# set up various models
occu_null <- occu(formula = ~ 1 ~1, data = cov_occu)
occu_m2 <- occu(formula = ~ temp + fasting_period ~1, data = cov_occu)
occu_m3 <- occu(formula = ~ temp * fasting_period ~1, data = cov_occu)
occu_m4 <- occu(formula = ~ fasting_period ~1, data = cov_occu)
occu_m5 <- occu(formula = ~ temp ~1, data = cov_occu)
occu_m6 <- occu(formula = ~ temp + fasting_period + weight_gain ~1, data = cov_occu)
occu_m7 <- occu(formula = ~ weight_gain ~1, data = cov_occu)
occu_m8 <- occu(formula = ~ temp + fasting_period + weight_gain + days_attached ~1, data = cov_occu)
occu_m9 <- occu(formula = ~ fasting_period + days_attached ~1, data = cov_occu)

# checking fit
fit <- fitList('psi(.)p(.)' = occu_null,
               'psi(.)p(temp + fasting_period)' = occu_m2,
               'psi(.)p(temp + fasting_period + temp*fasting_period)' = occu_m3,
               'psi(.)p(fasting_period)' = occu_m4,
               'psi(.)p(temp)' = occu_m5,
               'psi(.)p(temp + fasting_period + weight_gain)' = occu_m6,
               'psi(.)p(weight_gain)' = occu_m7,
               'psi(.)p(temp + fasting_period + weight_gain + days_attached)' = occu_m8,
               'psi(.)p(fasting_period + days_attached)' = occu_m9)
modSel(fit)

# looking at just fasting period model
preds1 <- predict(occu_m4, type ="det", new = data.frame(fasting_period = c(0:30)))
preds1$fasting_period <- seq(0, 30, by = 1)

ggplot(data = preds1, aes(x = fasting_period, y = Predicted)) +
  geom_smooth(stat = "smooth") +
  scale_y_continuous(limits = c(0.4, 1)) +
  geom_ribbon(data = preds1, aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(title = "Predicted Host Detection Probability",
       x = "Fasting Period (Days)",
       y = "Detection Probability") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold"), 
        axis.title.y = element_text(face = "bold"))




# looking at fasting period + days attached model

# back transform
backTransform(occu_m9, type = "state")

# generate a prediction grid for both variables
days_attached <- seq(0, 10, length.out = 100)
fasting_period <- seq(0, 30, length.out = 100)

# create an expanded data frame
pred_grid <- expand.grid(days_attached = days_attached, fasting_period = fasting_period)

# set up predictions
preds2 <- predict(occu_m9, type = "det", newdata = pred_grid)
pred_grid$preds <- preds2$Predicted # add to grid
# visualize data

# heatmap
ggplot(pred_grid, aes(x = fasting_period, y = days_attached, fill = preds)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "yellow", "red")) +
  labs(title = "Predicted Host Detection Probability",
       x = "Fasting Period (Days)",
       y = "Days Attached",
       fill = "Probability") +
  theme_minimal()

# contour plot
ggplot(pred_grid, aes(x = fasting_period, y = days_attached, z = preds)) +
  geom_contour_filled(breaks = seq(0, 1, by = 0.1)) +
  geom_textcontour(linecolour = "#242424", textcolour = "#242424") + 
  labs(title = "Predicted Host Detection Probability",
       x = "Fasting Period (Days)",
       y = "Days Attached",
       fill = "Probability") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12))




# Occupancy Model for a Single Replicate
detections_single <- data.frame(rep1 = full_detection_data$detection_rep1)
# Create the occupancy frame object with the modified detection data
cov_count_single <- unmarkedFramePCount(y = detections_single, siteCovs = site_covs)
summary(cov_count_single)

# Set up the model using the single replicate
count_m4_single <- pcount(formula = ~ fasting_period ~1, data = cov_count_single)

# Checking fit with only the single replicate model
fit_single <- fitList('lambda(.)p(fasting_period)' = count_m4_single)

# Make predictions with the single replicate model
preds_single <- predict(count_m4_single, type = "det", newdata = data.frame(fasting_period = c(0:30)))
preds_single$fasting_period <- seq(0, 30, by = 1)

# Plot predictions
ggplot(data = preds_single, aes(x = fasting_period, y = Predicted)) +
  geom_smooth(stat = "smooth") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(title = "Predicted Host Detection Probability with Single Replicate",
       x = "Fasting Period (Days)",
       y = "Detection Probability") +
  theme_minimal()


