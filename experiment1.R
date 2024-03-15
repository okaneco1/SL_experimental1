# Experimental Trials for Temperature/Fasting Period 

# libraries
library(tidyverse)
library(readxl)

# import data
com_mat <- read_excel("submission1A_community_matrix.xls", range = "A1:AW506")
dissection_data <- read_excel("dissection data.xlsx", range = "A1:N97")

# adjustments to data
colnames(com_mat)[1] <- "sample"

# limiting to just samples for experimental 2022 samples
com_mat <- com_mat %>%
  filter(grepl("T22", sample))

# add padding to sample names and sort
com_mat$sample <- sub("T22_(\\d)\\.", "T22_0\\1\\.", com_mat$sample)
com_mat <- com_mat[order(com_mat[[1]], decreasing = FALSE), ]

# remove columns with no sequences
for (i in ncol(com_mat):2) { # looped in reverse to avoid missing columns during the loop
  if (sum(com_mat[,i]) < 1) {
    com_mat[,i] <- NULL
  }
}

# only removed two columns

# next step, just get a sum of every column and display to see which species
# to remove based on what is known to be contamination



