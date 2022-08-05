################################################################################
# Load libraries
################################################################################

library(tidyverse)
library(icesTAF)
library(cowplot)
theme_set(theme_cowplot())

################################################################################
# Define study/data paths
################################################################################

github_dir <- "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Data_Prep_QC/")
icesTAF::mkdir(plot_path)
data_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/"
rdata_path <- paste0(data_path, "Rdata/")

source(paste0(github_dir, "helper_functions/QC_functions.R"))

input_dataset_name = "UCLA"

################################################################################
# Movement in scanner
################################################################################

movement_data <- compile_movement_data(fd_path = paste0(data_path, "movementData/"),
                                       input_dataset_name = input_dataset_name,
                                       subject_csv = paste0(data_path, "participants.csv")) %>%
  mutate(diagnosis = stringr::str_to_sentence(diagnosis))


# Plot movement data by diagnosis
plot_FD_vs_diagnosis(movement_data = movement_data,
                     input_dataset_name = input_dataset_name)
ggsave(paste0(plot_path, 
              sprintf("Fractional_Displacement_by_Subject_%s.png", 
                      input_dataset_name)),
       width=4, height=4, units="in", dpi=300)

# Plot number of subjects retained at each FD threshold
plot_subjects_per_fd_threshold(movement_data,
                               input_dataset_name = input_dataset_name)

ggsave(paste0(plot_path, 
              sprintf("Subjects_Retained_by_FD_Threshold_%s.png", 
                      input_dataset_name)),
       width=5, height=6, units="in", dpi=300)

# Plot Schz:Control ratio retained at each FD threshold
plot_schz_ctrl_ratio_per_fd_threshold(movement_data,
                                      input_dataset_name = input_dataset_name)

ggsave(paste0(plot_path, 
              sprintf("SCZ_to_CTRL_Retained_by_FD_Threshold_%s.png", 
                      input_dataset_name)),
       width=5, height=4, units="in", dpi=300)

################################################################################
# Univariate data
################################################################################

#### catch22 + catchaMouse16

#### Subjects with all NA values by univariate feature set
for (feature_set in c("catch22", "catchaMouse16")) {
  NA_subjects <- find_univariate_subject_na(rdata_path = rdata_path,
                                            feature_set = feature_set)
  
  plot_NA_subject_ts(rdata_path=rdata_path, 
                     feature_set = feature_set,
                     input_dataset_name = input_dataset_name,
                     NA_subject_IDs=NA_subjects$Subject_ID)
  
  ggsave(paste0(plot_path, 
                sprintf("Raw_fMRI_Signal_NA_Subjects_Univariate_%s.png", 
                        feature_set)),
         width=7, height=6, units="in", dpi=300)
}

#### Features that are NA for all subjects
for (feature_set in c("catch22", "catchaMouse16")) {
  NA_features <- find_univariate_feature_na(rdata_path = rdata_path,
                                            feature_set = feature_set)
  if (nrow(NA_features) > 0) {
    cat("\nNA features for", feature_set, "= ", unique(NA_features$names))
  }
}

# catchaMouse16 ST_LocalExtrema_n100_diffmaxabsmin is NA because it requires > 200 time points,
# Which the UCLA schizophrenia dataset does not have