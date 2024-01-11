library(tidyverse)
library(readr)
library(ggrepel)
library(ggstatsplot)
library(ggpubr)


# My data

my_dat <- read.csv(file = "../processed_data/NAguideR_seqknn.csv")
my_dat_mat <- my_dat[,-1] %>% as.matrix()
row.names(my_dat_mat) <- my_dat$X

# Validation data

val_dat <- read.csv("../processed_data/validation_data_seqknn.csv") 
val_dat <- separate_longer_delim(val_dat,X,";")
val_dat_mat <- val_dat[,-1] %>% as.matrix()
row.names(val_dat_mat) <- val_dat$X

# 
dim(my_dat_mat)
dim(val_dat_mat)
sum(rownames(my_dat_mat) %in% rownames(val_dat_mat))

com_proteins <- base::intersect(rownames(my_dat_mat),rownames(val_dat_mat))

my_dat_mat_all <- my_dat_mat[com_proteins,] %>% t()
val_dat_mat_all <- val_dat_mat[com_proteins,] %>% t()

## scale

my_scale <- 1/mean(my_dat_mat_all)
my_dat_mat_all <- my_dat_mat_all * my_scale

val_scale <- 1/mean(val_dat_mat_all)
val_dat_mat_all <- val_dat_mat_all * val_scale


# select difference proteins

my_diff_res <- read.csv("../results/diff_result.csv")
val_diff_res <- read.csv("../results/diff_val_result.csv")

m_dat <- merge(my_diff_res,val_diff_res,by = "Protein")

com_diff_proteins <- m_dat$Protein[m_dat$pvalue.x< 0.05 & m_dat$pvalue.y < 0.05 & (log2(m_dat$ratio.x) * log2(m_dat$ratio.y)) >0]

my_dat_mat_diff <- my_dat_mat_all[,com_diff_proteins]
val_dat_mat_diff <- val_dat_mat_all[,com_diff_proteins]

write.csv(my_dat_mat_diff,file = "../processed_data/my_dat_ml_diff.csv",row.names = F)
write.csv(val_dat_mat_diff,file = "../processed_data/val_dat_ml_diff.csv",row.names = F)

