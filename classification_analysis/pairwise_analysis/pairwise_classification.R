#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------
# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--github_dir", default="/headnode1/abry4213/github/")
parser$add_argument("--data_path", default="/headnode1/abry4213/data/UCLA_Schizophrenia/")
parser$add_argument("--sample_metadata_file", default="UCLA_Schizophrenia_sample_metadata.Rds")
parser$add_argument("--pairwise_feature_set", default="pyspi14")
parser$add_argument("--univariate_feature_set", default="catch22")
parser$add_argument("--input_feature_data_file")
parser$add_argument("--dataset_ID", default="UCLA_Schizophrenia")
parser$add_argument("--email")

# Parse input arguments
args <- parser$parse_args()
github_dir <- args$github_dir
data_path <- args$data_path
rdata_path <- args$rdata_path
pairwise_feature_set <- args$pairwise_feature_set
univariate_feature_set <- args$univariate_feature_set
input_feature_data_file <- args$input_feature_data_file
dataset_ID <- args$dataset_ID
sample_metadata_file <- args$sample_metadata_file
email <- args$email

# univariate_feature_set <- "catch22"
# pairwise_feature_set <- "pyspi14"
# github_dir <- "~/github/"
# email <- "abry4213@uni.sydney.edu.au"
# data_path <- "~/data/UCLA_CNP_ABIDE_ASD/"
# input_feature_data_file <- "UCLA_CNP_AROMA_2P_GMR_ABIDE_ASD_FC1000_pyspi14_filtered_zscored.Rds"
# dataset_ID <- "UCLA_CNP_ABIDE_ASD"
# sample_metadata_file <- "UCLA_CNP_ABIDE_ASD_sample_metadata.Rds"

rdata_path <- paste0(data_path, "processed_data/Rdata/")
plot_dir <- paste0(data_path, "plots/")

TAF::mkdir(plot_dir)

# Set the seed
set.seed(127)

# Load tidyverse
library(tidyverse)

#-------------------------------------------------------------------------------
# Source helper scripts
#-------------------------------------------------------------------------------
helper_script_dir = paste0(github_dir, "fMRI_FeaturesDisorders/helper_functions/classification/")
source(paste0(helper_script_dir, "Linear_SVM.R"))
source(paste0(helper_script_dir, "Null_distributions.R"))

subjects_to_use <- readRDS(paste0(rdata_path, sprintf("%s_samples_with_univariate_%s_and_pairwise_%s_filtered.Rds",
                                                      dataset_ID,
                                                      univariate_feature_set,
                                                      pairwise_feature_set)))

# Load sample metadata
sample_metadata <- readRDS(paste0(data_path, "study_metadata/", sample_metadata_file)) %>%
  filter(Sample_ID %in% subjects_to_use)

################################################################################
# Load pyspi data
################################################################################
pairwise_feature_data <- readRDS(paste0(rdata_path, input_feature_data_file)) %>%
  mutate(Diagnosis = ifelse(Diagnosis == "BIPOLAR", "Bipolar", Diagnosis))

SPI_directionality_file <- paste0(github_dir, "fMRI_FeaturesDisorders/classification_analysis/pairwise_analysis/SPI_Direction_Info.csv")
SPI_directionality <- read.csv(SPI_directionality_file)

# Load case data
pairwise_features_groups <- pairwise_feature_data %>%
  left_join(., sample_metadata) %>%
  filter(Diagnosis != "Control")

# Load control data
pairwise_features_controls <- pairwise_feature_data %>%
  left_join(., sample_metadata) %>%
  filter(Diagnosis == "Control")


################################################################################
# Create ten folds to use for all analyses
################################################################################
# Load contrasts to control participants
study_comparisons_to_control <- readRDS(paste0(rdata_path, dataset_ID, "_comparisons_to_control.Rds"))

# Use 10-fold CV
k <- 10

# Iterate over each contrast to controls per dataset
if (!file.exists(paste0(rdata_path, dataset_ID, "_samples_per_10_folds_10_repeats.Rds"))) {
  sample_folds_list <- list()
  for (i in 1:nrow(study_comparisons_to_control)) {
    study <- study_comparisons_to_control$Study[i]
    group_to_compare <- study_comparisons_to_control$Group_to_Compare[i]
    
    group_folds <- list()
    
    control_data <- sample_metadata %>%
      filter(Study == study, 
             Diagnosis == "Control") %>%
      distinct(Sample_ID, Diagnosis) 
    
    control_vector <- control_data %>%
      pull(Diagnosis)
    
    group_data <- sample_metadata %>%
      filter(Study == study,
             Diagnosis == group_to_compare) %>%
      distinct(Sample_ID, Diagnosis) 
    
    group_vector <- group_data %>%
      pull(Diagnosis)
    
    merged_data <- plyr::rbind.fill(control_data,
                                    group_data)
    
    # 10 repeats
    for (r in 1:10) {
      # 10 folds per repeat for 10-fold CV
      group_folds_r <- caret::createFolds(c(control_vector, group_vector),
                                          k = k,
                                          list = TRUE,
                                          returnTrain = FALSE)
      
      # Map indices to sample IDs to avoid ambiguity
      group_folds_samples_r <- purrr::map(group_folds_r, ~ (merged_data[.x,]) %>% pull(Sample_ID))
      
      # Append to list
      group_folds[[r]] <- group_folds_samples_r
    }
    
    # Combine group splits
    group_df <- as.data.frame(cbind(Study = rep(study, 10),
                                    Group_to_Compare = rep(group_to_compare, 10),
                                    Repeat_Number = 1:10,
                                    Folds = group_folds))
    
    # Append to dataframe
    sample_folds_list <- list.append(sample_folds_list, group_df)
  }
  
  # Combine fold results into a dataframe
  sample_folds <- do.call(plyr::rbind.fill, sample_folds_list)
  
  # Save to an Rds file
  saveRDS(sample_folds, paste0(rdata_path, dataset_ID, 
                               sprintf("_samples_per_%s_folds_10_repeats.Rds",
                                       k)))
} else {
  sample_folds <- readRDS(paste0(rdata_path, dataset_ID, 
                                 sprintf("_samples_per_%s_folds_10_repeats.Rds",
                                         k)))
}


################################################################################
# Define weighting parameters
################################################################################
# Use a linear kernel
svm_kernel = "linear"
weighting_name <- "inv_prob"
use_inv_prob_weighting <- TRUE

grouping_param_df <- data.frame(grouping_var = c("SPI"),
                                SVM_feature_var = c("region_pair"))

output_file_RDS <- paste0(rdata_path, 
                          sprintf("%s_pairwise_CV_linear_SVM_%s_%s.Rds",
                                  dataset_ID,
                                  pairwise_feature_set, 
                                  weighting_name))

if (!file.exists(output_file_RDS)) {
  # Load non-control data 
  
  pairwise_class_res_list <- list()
  
  for (i in 1:nrow(study_comparisons_to_control)) {
    study <- study_comparisons_to_control$Study[i]
    group_to_compare <- study_comparisons_to_control$Group_to_Compare[i] 
    
    # Get subjects to compare
    sample_groups <- sample_metadata %>%
      filter(Sample_ID %in% subjects_to_use,
             Study == study,
             Diagnosis %in% c("Control", group_to_compare)) %>%
      distinct(Sample_ID, Diagnosis)
    
    # Isolate group folds
    group_folds <- sample_folds %>%
      filter(Study == study,
             Group_to_Compare == group_to_compare) %>%
      pull(Folds)
    
    # Combine control data with diagnosis group data
    group_data_for_SVM <- pairwise_features_groups %>%
      filter(Study == study,
             Diagnosis == group_to_compare) %>%
      plyr::rbind.fill(., pairwise_features_controls %>%
                         filter(Study == study)) 
    
    group_data_for_SVM_feature_set <- group_data_for_SVM %>%
      filter(feature_set == pairwise_feature_set)
    
    # Run SVM
    for (j in 1:nrow(grouping_param_df)) {
      grouping_var = grouping_param_df$grouping_var[j]
      SVM_feature_var = grouping_param_df$SVM_feature_var[j]
      
      group_wise_SVM_CV_feature_set <- 1:length(group_folds) %>%
        purrr::map_df( ~ run_pairwise_cv_svm_by_input_var(feature_matrix = group_data_for_SVM_feature_set,
                                                          dataset_ID = dataset_ID,
                                                          pairwise_feature_set = pairwise_feature_set,
                                                          sample_groups = sample_groups,
                                                          SPI_directionality = SPI_directionality,
                                                          svm_kernel = svm_kernel,
                                                          grouping_var = grouping_var,
                                                          svm_feature_var = SVM_feature_var,
                                                          flds = group_folds[[.x]],
                                                          repeat_number = .x,
                                                          out_of_sample_only = TRUE,
                                                          drop_NaN = FALSE,
                                                          impute_NaN = TRUE,
                                                          use_inv_prob_weighting = use_inv_prob_weighting)) %>%
        mutate(Study = study,
               Analysis_Type = grouping_type,
               Group_to_Compare = group_to_compare)
      
      pairwise_class_res_list <- list.append(pairwise_class_res_list, group_wise_SVM_CV_feature_set)
    }
  }
  pairwise_class_res <- do.call(plyr::rbind.fill, pairwise_class_res_list)
  saveRDS(pairwise_class_res, file = output_file_RDS)
} else {
  pairwise_class_res <- readRDS(file = output_file_RDS)
}

