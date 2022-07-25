github_dir <- "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/"
data_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/"
rdata_path <- paste0(data_path, "Rdata/")
set.seed(127)

library(tidyverse)

# Load info on univariate (theft) subjects
input_dataset_name = "UCLA"
feature_set = "catch22"
univariate_subject_info <- readRDS(paste0(rdata_path, "UCLA_filtered_subject_info_catch22.Rds"))
pairwise_subject_info <- readRDS(paste0(rdata_path, "Filtered_subject_info_pyspi_19.Rds"))

intersection <- inner_join(pairwise_subject_info, univariate_subject_info)
saveRDS(intersection, file=paste0(rdata_path, "UCLA_Subjects_with_Univariate_and_Pairwise.Rds"))
