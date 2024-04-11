# Does Variance Increase as a Function of Temp and Fasting Period?

# libraries
library(tidyverse)
library(readxl)

#----- community data
com_mat <- read_excel("replicate_calculations/submission1A_community_matrix.xls", range = "A1:AW506")
com_mat_r <- read_excel("replicate_calculations/submission1B_community_matrix.xls")


#----- organizing data

colnames(com_mat)[1] <- "sample"
colnames(com_mat_r)[1] <- "sample"

# making all names the same
com_mat$sample <- sub("\\.12S", "", com_mat$sample)
com_mat_r$sample <- sub("_R\\.12S", "", com_mat_r$sample)

# reorder
com_mat <- com_mat[order(com_mat[[1]], decreasing = FALSE), ]
com_mat_r <- com_mat_r[order(com_mat_r[[1]], decreasing = FALSE), ]

# limit to just experimental samples
com_mat <- filter(com_mat, grepl("T2", com_mat$sample))
com_mat_r <- filter(com_mat_r, grepl("T2", com_mat_r$sample))

# function to transform/unify the sample names (to match dissection data)
transform_sample_names <- function(name) {
  # remove the 12S
  name <- gsub("\\.\\d+S$", "", name)
  # add 0 before single digits
  name <- sub("[_-](\\d)$", "-0\\1", name)
  # replace _ with -
  name <- sub("_", "-", name)
  # insert 20 after T
  name <- sub("T", "T20", name)
  return(name)
}
com_mat$sample <- sapply(com_mat$sample, transform_sample_names)
com_mat_r$sample <- sapply(com_mat_r$sample, transform_sample_names)


#----- merged replicate data frames on the basis of the abs. value of the difference

merged_df <- full_join(com_mat, com_mat_r, by = "sample", suffix = c(".o", ".r"))


# initialize data frame
abs_val_df <- data.frame(sample = merged_df[,1])

for (i in 2:ncol(merged_df)) {
  # start with ".o" columns and find absolute difference between ".r" columns
  if (grepl(".o$", colnames(merged_df[i]))) {
    # get names, assign columns
    otu_name <- sub("\\.o", "", colnames(merged_df[i]))
    o_column <- merged_df[,i]
    r_column_name <- paste0(otu_name, ".r")
    
    # check if the replicate column exists
    if (!r_column_name %in% colnames(merged_df)) {
      # error if .r column is missing
      stop(paste("The corresponding .r column for", otu_name, "is missing."))
    }
    
    r_column <- merged_df[, r_column_name]
    
    # set NAs to 0, and put in data frame
    r_column[is.na(r_column)] <- 0
    o_column[is.na(o_column)] <- 0
    
    # calculate absolute value of differences
    combined_df <- data.frame(o = o_column, r = r_column)
    abs_diff <- pull(abs(o_column - r_column))
    
    # add to data frame
    abs_val_df <- abs_val_df %>%
      mutate("{otu_name}_abs_diff" := abs_diff)
  } else {
    # if column does not end in ".r" either (only reads from one batch)
    if (!grepl(".r$", colnames(merged_df[i]))) {
      otu_name <- colnames(merged_df[i])
      
      # abs. diff. will be equal to the read count, as no replicate exists
      abs_diff <- pull(merged_df[,i])
      abs_diff[is.na(abs_diff)] <- 0 # NAs to 0
      
      # add to data frame
      abs_val_df <- abs_val_df %>%
        mutate("{otu_name}_abs_diff" := abs_diff)
    } else {
      # if neither of first options, column must end in ".r", so skip it 
      # as it is already accounted for
      next
    }
  }
}

view(abs_val_df)




#------- trying with log transformation
#----- merged replicate data frames on the basis of the abs. value of the difference

merged_df <- full_join(com_mat, com_mat_r, by = "sample", suffix = c(".o", ".r"))


# initialize data frame
log_abs_val_df <- data.frame(sample = merged_df[,1])

for (i in 2:ncol(merged_df)) {
  # start with ".o" columns and find absolute difference between ".r" columns
  if (grepl(".o$", colnames(merged_df[i]))) {
    # get names, assign columns
    otu_name <- sub("\\.o", "", colnames(merged_df[i]))
    o_column <- merged_df[,i]
    r_column_name <- paste0(otu_name, ".r")
    
    # check if the replicate column exists
    if (!r_column_name %in% colnames(merged_df)) {
      # error if .r column is missing
      stop(paste("The corresponding .r column for", otu_name, "is missing."))
    }
    
    r_column <- merged_df[, r_column_name]
    
    # set NAs to 0, and put in data frame
    r_column[is.na(r_column)] <- 0
    o_column[is.na(o_column)] <- 0
    
    # calculate absolute value of differences
    # combined_df <- data.frame(log_o = log(o_column), log_r = log(r_column))
    
    # take log of columns
    o_log <- log(o_column)
    r_log <- log(r_column)
    
    # change -Inf values to 0
    o_log <- o_log %>%
      mutate(across(1, ~if_else(. == -Inf, 0, .)))
    r_log <- r_log %>%
      mutate(across(1, ~if_else(. == -Inf, 0, .)))
    
    # take absolute difference of logs
    log_abs_diff <- pull(abs(o_log - r_log))
    
    # add to data frame
    log_abs_val_df <- log_abs_val_df %>%
      mutate("{otu_name}_log_abs_diff" := log_abs_diff)
  } else {
    # if column does not end in ".r" either (only reads from one batch)
    if (!grepl(".r$", colnames(merged_df[i]))) {
      otu_name <- colnames(merged_df[i])
      
      # log abs. diff. will be equal to the log of the read count, as no replicate exists
      log_abs_diff <- pull(log(merged_df[,i]))
      log_abs_diff[is.na(abs_diff)] <- 0 # NAs to 0
      
      # add to data frame
      log_abs_val_df <- log_abs_val_df %>%
        mutate("{otu_name}_log_abs_diff" := log_abs_diff)
    } else {
      # if neither of first options, column must end in ".r", so skip it 
      # as it is already accounted for
      next
    }
  }
}

view(log_abs_val_df)





#----- add in dissection data and merge
# read in
dissection_data_5C <- read_excel("Sea Lamprey Feeding Study(2022_ROB_441008)_TubeID_CO.xlsx", range = "A1:L27", sheet = 1)
dissection_data_10C <- read_excel("Sea Lamprey Feeding Study(2022_ROB_441008)_TubeID_CO.xlsx", range = "A1:L28", sheet = 2)
dissection_data_15C <- read_excel("Sea Lamprey Feeding Study(2022_ROB_441008)_TubeID_CO.xlsx", range = "A1:L27", sheet = 3)

# change "Tube ID" column name to "sample" for merge
colnames(dissection_data_5C)[1] <- "sample"
colnames(dissection_data_10C)[1] <- "sample"
colnames(dissection_data_15C)[1] <- "sample"

# add temp column and bind
dissection_data_5C <- mutate(dissection_data_5C, temp = 5, .after = 4)
dissection_data_10C <- mutate(dissection_data_10C, temp = 10, .after = 4)
dissection_data_15C <- mutate(dissection_data_15C, temp = 15, .after = 4)
dissection_data <- rbind(dissection_data_5C, dissection_data_10C, dissection_data_15C)

# merge with community matrix
full_av_data <- inner_join(abs_val_df, dissection_data, by = "sample")
full_log_av_data <- inner_join(log_abs_val_df, dissection_data, by = "sample")

# data frame adjustments (abs diff)
colnames(full_av_data)[63] <- "fasting_period"
full_av_data <- full_av_data[ ,c(1:9, 54:ncol(full_av_data))] # removing low count OTUs
full_av_data <- as.data.frame(full_av_data)
# all trout
full_av_data <- full_av_data %>%
  mutate(all_trout_ad = Salvelinus_namaycush_abs_diff + Salmonidae_unclassified_abs_diff)

# data frame adjustments (log abs diff)
colnames(full_log_av_data)[63] <- "fasting_period"
full_log_av_data <- full_log_av_data[ ,c(1:9, 54:ncol(full_log_av_data))] # removing low count OTUs
# all trout
full_log_av_data <- full_log_av_data %>%
  mutate(all_trout_log_ad = Salvelinus_namaycush_log_abs_diff + Salmonidae_unclassified_log_abs_diff)





#-------- comparisons

# viewing data
sorted_data <- full_av_data %>%
  dplyr::select(sample, Salvelinus_namaycush_abs_diff, temp, fasting_period) %>%
  dplyr::arrange(temp, fasting_period)


# visualizing
ggplot(full_log_av_data, aes(x = fasting_period, y = all_trout_log_ad, color = factor(temp))) +
  geom_boxplot(aes(group = interaction(fasting_period, factor(temp))),
               alpha = 0.5, position = position_dodge(width = 3), 
               outlier.shape = NA,  coef = Inf) +
  geom_smooth(method = "lm", aes(group = factor(temp)), formula = y ~ x, se = FALSE) +
  labs(color = "Temperature (°C)", y = "Log Absolute Dif. Between Lake Trout Replicates", x = "Fasting Period (Days)") +
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_brewer(palette = "Dark2") +
  scale_x_continuous(breaks = c(0, 5, 10, 20, 30)) +
  facet_wrap(vars(temp), labeller = labeller(temp = custom_labels)) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold", size = 12))



#------- statistical comparisons

#--- ANOVAS
aov <- aov(all_trout_log_ad ~ as.factor(temp) * as.factor(fasting_period), data = full_log_av_data) 

summary(aov)

#---- testing assumptions
library(car)
#verifying normality
plot(aov, which = 2)
qqPlot(aov$residuals, id = FALSE)

hist(resid(aov), breaks = 10)

# homogeneity of variances
plot(aov, which = 3)







