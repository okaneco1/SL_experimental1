# Combining Both Replicate and Original Datasets

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

# add padding for zeros (DON'T REALLY NEED THIS)
#com_mat$sample <- sub("_(\\d)$", "_0\\1", com_mat$sample)
#com_mat$sample <- sub("_([A-Z])(\\d)$", "_\\10\\2", com_mat$sample)
#com_mat_r$sample <- sub("_(\\d)$", "_0\\1", com_mat_r$sample)
#com_mat_r$sample <- sub("_[A-Z](\\d)$", "_0\\1", com_mat_r$sample)

# reorder
com_mat <- com_mat[order(com_mat[[1]], decreasing = FALSE), ]
com_mat_r <- com_mat_r[order(com_mat_r[[1]], decreasing = FALSE), ]



#----- combining data
merged_df <- full_join(com_mat, com_mat_r, by = "sample", suffix = c(".o", ".r"))


# initialize data frame
average_df <- data.frame(sample = merged_df[,1])

for (i in 2:ncol(merged_df)) {
  # start with ".o" columns and average with ".r" columns
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
    
    # calculate means
    combined_df <- data.frame(o = o_column, r = r_column)
    mean_column <- rowMeans(combined_df)
    
    # add to data frame
    average_df <- average_df %>%
      mutate("{otu_name}_mean" := mean_column)
  } else {
    # get columns that do not end in ".r" either (only reads from one batch)
    if (!grepl(".r$", colnames(merged_df[i]))) {
      otu_name <- colnames(merged_df[i])
      
      # average values with 0 (as replicate was non-existent, thus 0 reads)
      mean_column <- merged_df[,i] / 2
      mean_column[is.na(mean_column)] <- 0 # NAs to 0
      
      # add to data frame
      average_df <- average_df %>%
        mutate("{otu_name}_mean" := mean_column[,1])
    } else {
      # column must end in ".r", so skip it
      next
    }
  }
}

view(average_df)


#---- write out
write.csv(average_df, "replicate_calculations/submission1_averaged_community_matrix.csv", row.names = FALSE)


#---- which samples were only present in one of the datasets?
sample_test <- data.frame(sample = c(com_mat$sample, com_mat_r$sample))
unique_samples <- sample_test %>%
  group_by(sample) %>%
  summarise(Count = n()) %>%
  filter(Count == 1) %>%
  select(sample)

# To see the results
print(unique_samples)



