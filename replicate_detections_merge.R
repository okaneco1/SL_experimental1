# Categorize Detections Between Replicates

# The objective here is to determine which replicates detected lake trout
# sequences, and to compare whether double, single, or no detections
# correlates with differences in temperature and fasting period

#----- libraries
library(tidyverse)
library(readxl)
library(patchwork)
library(sjPlot)
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


# set categorical variables as factors 
full_detection_data$temp <- factor(full_detection_data$temp,
                                   levels = c(5, 10, 15))
full_detection_data$fasting_period <- factor(full_detection_data$fasting_period,
                                             levels = c(0, 5, 10, 20, 30))
full_detection_data$combined_detection <- factor(full_detection_data$combined_detection,
                                                 levels = c(0, 1, 2))

# create additional data column for both detections (and label)
full_detection_data$both_det <- ifelse(full_detection_data$combined_detection == 2, 1, 0)
full_detection_data$both_det <- factor(full_detection_data$both_det,
                                       levels = c(0, 1),
                                       labels = c("no detection", "detection"))


#------ can now look for patterns where detection is in *at least one* replicate
# create new columns (and label)
full_detection_data$one_det <- ifelse(full_detection_data$combined_detection == 1 |
                                        full_detection_data$combined_detection == 2, 1, 0)
full_detection_data$one_det <- factor(full_detection_data$one_det,
                                      levels = c(0, 1),
                                      labels = c("no detection", "detection"))



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


#--------------------------------------------------------
# OCCUPANCY MODEL
#--------------------------------------------------------

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
occu_m1 <- occu(formula = ~ temp ~ 1, data = cov_occu)
occu_m2 <- occu(formula = ~ fasting_period ~ 1, data = cov_occu)
occu_m3 <- occu(formula = ~ days_attached ~ 1, data = cov_occu)
occu_m4 <- occu(formula = ~ weight_gain ~ 1, data = cov_occu)
occu_m5 <- occu(formula = ~ temp + fasting_period ~ 1, data = cov_occu)
occu_m6 <- occu(formula = ~ temp + days_attached ~ 1, data = cov_occu)
occu_m7 <- occu(formula = ~ temp + weight_gain ~ 1, data = cov_occu)
occu_m8 <- occu(formula = ~ fasting_period + days_attached ~ 1, data = cov_occu)
occu_m9 <- occu(formula = ~ fasting_period + weight_gain ~ 1, data = cov_occu)
occu_m10 <- occu(formula = ~ days_attached + weight_gain ~ 1, data = cov_occu)
occu_m11 <- occu(formula = ~ temp + fasting_period + days_attached ~ 1, data = cov_occu)
occu_m12 <- occu(formula = ~ temp + fasting_period + weight_gain ~ 1, data = cov_occu)
occu_m13 <- occu(formula = ~ temp + days_attached + weight_gain ~ 1, data = cov_occu)
occu_m14 <- occu(formula = ~ fasting_period + days_attached + weight_gain ~ 1, data = cov_occu)
occu_m15 <- occu(formula = ~ temp + fasting_period + days_attached + weight_gain ~ 1, data = cov_occu)
occu_m16 <- occu(formula = ~ temp * fasting_period ~ 1, data = cov_occu)

# checking fit
fit <- fitList(
  'psi(.)p(.)' = occu_null,  # null model
  
  # single variable models
  'psi(.)p(temp)' = occu_m1,
  'psi(.)p(fasting_period)' = occu_m2,
  'psi(.)p(days_attached)' = occu_m3,
  'psi(.)p(weight_gain)' = occu_m4,
  # two variable models
  'psi(.)p(temp + fasting_period)' = occu_m5,
  'psi(.)p(temp + days_attached)' = occu_m6,
  'psi(.)p(temp + weight_gain)' = occu_m7,
  'psi(.)p(fasting_period + days_attached)' = occu_m8,
  'psi(.)p(fasting_period + weight_gain)' = occu_m9,
  'psi(.)p(days_attached + weight_gain)' = occu_m10,
  # three variable models
  'psi(.)p(temp + fasting_period + days_attached)' = occu_m11,
  'psi(.)p(temp + fasting_period + weight_gain)' = occu_m12,
  'psi(.)p(temp + days_attached + weight_gain)' = occu_m13,
  'psi(.)p(fasting_period + days_attached + weight_gain)' = occu_m14,
  # four variable model
  'psi(.)p(temp + fasting_period + days_attached + weight_gain)' = occu_m15,
  # interaction model
  'psi(.)p(temp * fasting_period)' = occu_m16
)
mod_sel <- modSel(fit)

# write out model selection
modSel_df <- as.data.frame(mod_sel@Full)
write.csv(modSel_df, file = "model selection table occupancy exp1.csv", row.names = FALSE)

# looking at just fasting period model
preds1 <- predict(occu_m2, type ="det", new = data.frame(fasting_period = c(0:30)))
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


1 - exp(-0.0395) # fasting days
exp(0.5362) # days attached

# looking at fasting period + days attached model

# back transform
backTransform(occu_m8, type = "state")

# generate a prediction grid for both variables
days_attached <- seq(0, 10, length.out = 100)
fasting_period <- seq(0, 30, length.out = 100)

# create an expanded data frame
pred_grid <- expand.grid(days_attached = days_attached, fasting_period = fasting_period)

# set up predictions
preds2 <- predict(occu_m8, type = "det", newdata = pred_grid)
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
  labs(x = "Fasting Period (Days)",
       y = "Days Attached",
       fill = "Predicted Host\nDNA Detection\nProbability") +
  scale_x_continuous(limits = c(1, NA), breaks = seq(1, 30, by = 4), expand = c(0, 0)) +  # Start x-axis at 1
  scale_y_continuous(limits = c(NA, 7), expand = c(0, 0)) +  # End y-axis at 7
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12))

