# Comparing Community Matrix Output with Replicate Data

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

# lake trout specific plot
lt_df <- data.frame(sample = com_mat_r[,1], replicate_reads = com_mat_r$Salvelinus_namaycush)
lt_df <- lt_df %>% 
  left_join(com_mat[,c("sample", "Salvelinus_namaycush")], by = "sample")
lt_df[] <- lapply(lt_df, function(x) ifelse(is.na(x), 0, x)) # replace NAs with 0s

# plot
ggplot(lt_df, aes(x = Salvelinus_namaycush, y = replicate_reads))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0, color = "darkred")


#----- set up as a function

replicate_reads_plot <- function(column_name) {
  # ensure the column name is a character string
  column_name <- as.character(substitute(column_name))
  # set up data frame
  reads_df <- data.frame(sample = com_mat_r[,1], replicate_reads = com_mat_r[[column_name]])
  reads_df <- left_join(reads_df, com_mat[,c("sample", column_name)], by = "sample")
  # change NAs to 0s
  reads_df[] <- lapply(reads_df, function(x) ifelse(is.na(x), 0, x))
  # correlation test
  cor_test <- round(cor(reads_df[[column_name]], reads_df$replicate_reads, use = "complete.obs"), 3)
  # plotting
  p <- ggplot(reads_df, aes_string(x = column_name, y = "replicate_reads")) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "darkred") +
    geom_abline(slope = 2, intercept = 0, color = "darkgreen", linetype = "dashed", alpha = 0.5) +
    geom_abline(slope = 0.5, intercept = 0, color = "darkgreen", linetype = "dashed", alpha = 0.5) +
    labs(x = "Original Reads", y = "Replicate Reads", title = column_name)+
    theme_minimal() +
    annotate("text", x = Inf, y = Inf, 
             label = paste("R =", cor_test), 
             hjust = 2, vjust = 1, size = 5, color = "deepskyblue")
  return(p)
}



#----- testing with other OTUs

replicate_reads_plot(Salvelinus_namaycush)
replicate_reads_plot(Salmonidae_unclassified)
replicate_reads_plot(Catostomus_commersonii)
replicate_reads_plot(Petromyzontidae_unclassified)
replicate_reads_plot(Oncorhynchus_mykiss)
replicate_reads_plot(Lota_lota)
replicate_reads_plot(Catostomus_catostomus)

# combining plots
library(patchwork)

su <- replicate_reads_plot(Salmonidae_unclassified)
cc <- replicate_reads_plot(Catostomus_commersonii)
pu <- replicate_reads_plot(Petromyzontidae_unclassified)
om <- replicate_reads_plot(Oncorhynchus_mykiss)

combined_plot <- su + cc + pu + om
combined_plot





