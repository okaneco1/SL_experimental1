# Experimental Trials for Temperature/Fasting Period 

options(scipen = 999)

# libraries
library(tidyverse)
library(readxl)


# ---------- Data Imports and Cleaning

# import data
com_mat <- read_excel("submission1_averaged_community_matrix.xlsx", range = "A1:BA533")
dissection_data_5C <- read_excel("Sea Lamprey Feeding Study(2022_ROB_441008)_TubeID_CO.xlsx", range = "A1:L27", sheet = 1)
dissection_data_10C <- read_excel("Sea Lamprey Feeding Study(2022_ROB_441008)_TubeID_CO.xlsx", range = "A1:L28", sheet = 2)
dissection_data_15C <- read_excel("Sea Lamprey Feeding Study(2022_ROB_441008)_TubeID_CO.xlsx", range = "A1:L27", sheet = 3)

# adjustments to data
colnames(com_mat)[1] <- 'Tube ID'
com_mat <- filter(com_mat, grepl("T2", com_mat$`Tube ID`)) # limit to experimental samples

# align 'Tube ID' name structures
# function to transform the sample names
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

com_mat$`Tube ID` <- sapply(com_mat$`Tube ID`, transform_sample_names)

# add sample number and order
com_mat$sample <- gsub("T\\d{4}-", "", com_mat$`Tube ID`)
com_mat <- com_mat[, c(1, ncol(com_mat), 2:(ncol(com_mat)-1))] # places it as 2nd column
com_mat <- arrange(com_mat, sample)

# remove columns with no sequences
for (i in ncol(com_mat):3) { # looped in reverse to avoid missing columns during the loop
  if (sum(com_mat[,i]) < 1) {
    com_mat[,i] <- NULL
  }
}
# removed 19 columns


#---------- Merge Data
# add temperature to dissection data and combine
dissection_data_5C <- mutate(dissection_data_5C, temp = 5, .after = 4)
dissection_data_10C <- mutate(dissection_data_10C, temp = 10, .after = 4)
dissection_data_15C <- mutate(dissection_data_15C, temp = 15, .after = 4)
dissection_data <- rbind(dissection_data_5C, dissection_data_10C, dissection_data_15C)

# merge with community matrix
full_data <- inner_join(com_mat_filtered, dissection_data, by = 'Tube ID')
#full_data <- full_data[-(which(grepl("2023", full_data$`Tube ID`))), ] # remove single 2023 experimental sample

# write out this table for use in other scripts:
full_data <- as.data.frame(full_data)
write.csv(full_data, "Experiment_1_full_data_table.csv", row.names = F)


#---------- DATA VISUALIZATION
# make combined variable for temperature/fasting period
full_data <- full_data %>%
  mutate(temp_fast = paste(temp, `Fasting period(days)`, sep = "_"))

# order for plots
custom_order <- c("5_0", "5_5", "5_10", "5_20", "5_30",
                  "10_0", "10_5", "10_10", "10_20", "10_30",
                  "15_0", "15_5", "15_10", "15_20", "15_30")
#check
view(filter(full_data, grepl("^5_0$", full_data$temp_fast)))

# removing "_mean" from column names
colnames(full_data) <- sub("_mean", "", colnames(full_data))

#--- single variable comparisons

ggplot(full_data, aes(x = temp, y = Salvelinus_namaycush, group = `Fasting period(days)`)) +
  geom_boxplot()
ggplot(full_data, aes(x = temp, y = Salvelinus_namaycush, group = temp)) +
  geom_boxplot()


#--- both temperature and fasting period

# plot with lake trout only reads 
lake_trout_boxplot <- ggplot(full_data, aes(x = temp_fast, y = Salvelinus_namaycush)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  #geom_text(aes(label = `Tube ID`), position = position_jitter(width = 0.2, height = 0), hjust = 0, vjust = 0, check_overlap = TRUE) +
  scale_x_discrete(limits = custom_order)
lake_trout_boxplot

# allowing salmonidae unclassified as well
full_data <- full_data %>%
  mutate(trout_sum = Salvelinus_namaycush + Salmonidae_unclassified)
# plot
trout_box <- ggplot(full_data, aes(x = temp_fast, y = trout_sum)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_x_discrete(limits = custom_order)
trout_box







