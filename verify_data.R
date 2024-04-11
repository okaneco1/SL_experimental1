# Verify Temp and Fasting Period for each Tube

#--- libraries

library(tidyverse)
library(readxl)

#--- data

# dissection data
dissection_data <- read_excel("dissection data_highlighted.xlsx")

# feeding data
feeding_data_5C <- read_excel("Sea Lamprey Feeding Study(2022_ROB_441008)_TubeID_CO.xlsx", range = "A1:L27", sheet = 1)
feeding_data_10C <- read_excel("Sea Lamprey Feeding Study(2022_ROB_441008)_TubeID_CO.xlsx", range = "A1:L28", sheet = 2)
feeding_data_15C <- read_excel("Sea Lamprey Feeding Study(2022_ROB_441008)_TubeID_CO.xlsx", range = "A1:L27", sheet = 3)

# limit to experimental samples
com_mat <- filter(com_mat, grepl("T2", com_mat$`Tube ID`)) 

# merge feeding data
feeding_data_5C <- mutate(feeding_data_5C, temp = 5, .after = 4)
feeding_data_10C <- mutate(feeding_data_10C, temp = 10, .after = 4)
feeding_data_15C <- mutate(feeding_data_15C, temp = 15, .after = 4)
feeding_data <- rbind(feeding_data_5C, feeding_data_10C, feeding_data_15C)

# name Tube Id to sample
colnames(feeding_data)[1] <- "sample"
colnames(dissection_data)[1] <- "sample"

# put in order
feeding_data <- arrange(feeding_data, sample)
dissection_data <- arrange(dissection_data, sample)

# change names in dissection data to match feeding data
colnames(dissection_data)[8] <- "temp"
colnames(feeding_data)[11] <- "fasting"
colnames(dissection_data)[13] <- "fasting"


#---- verify if temp and fasting period in feeding data match dissection data
# set up data frame to capture results
verify_results <- data.frame(sample = feeding_data$sample, temp_match = NA, fast_match = NA)

# check 
for (i in 1:nrow(dissection_data)) {
  if (dissection_data$sample[i] %in% feeding_data$sample) {
    # get row for feeding data
    matching_row <- which(feeding_data$sample == dissection_data$sample[i])
    # does temp match?
    temp_match <- feeding_data$temp[matching_row] == dissection_data$temp[i]
    # add results to table
    verify_results[matching_row, "temp_match"] <- temp_match
    # does fasting period match?
    fast_match <- feeding_data$fasting[matching_row] == dissection_data$fasting[i]
    # add results to table
    verify_results[matching_row, "fast_match"] <- fast_match
  } else {
    # skip if sample is not in both datasets
    next
  }
}
verify_results

# two instances of false [fixed]





