# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--project_path", default="/project/hctsa/annie/")
parser$add_argument("--github_dir", default="/project/hctsa/annie/github/fMRI_FeaturesDisorders/")
parser$add_argument("--rdata_path", default="/project/hctsa/annie/data/scz/UCLA/Rdata/")
parser$add_argument("--pydata_path", default="/project/hctsa/annie/data/scz/UCLA/pydata/")
parser$add_argument("--feature_set", default="pyspi_19")
# project_path <- "D:/Virtual_Machines/Shared_Folder/"
# github_dir <- "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/"
# rdata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/Rdata/"
# pydata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/pydata/"
# feature_set <- "pyspi_19"

# Parse input arguments
args <- parser$parse_args()
project_path <- args$project_path
github_dir <- args$github_dir
rdata_path <- args$rdata_path
pydata_path <- args$pydata_path
feature_set <- args$feature_set

### Source functions
# Main
source(paste0(github_dir, "helper_functions/Linear_SVM.R"))
source(paste0(github_dir, "helper_functions/Visualization.R"))
source(paste0(github_dir, "helper_functions/Null_distributions.R"))

set.seed(127)

# Compare just AROMA+2P+GMR
noise_proc = "AROMA+2P+GMR"

# Use e1071 SVM with a linear kernel
test_package = "e1071"
kernel = "linear"

# ###############################################################################
# Load pyspi data
# ###############################################################################
pyspi_data_file <- paste0(pydata_path, "UCLA_all_subject_pyspi_19_AROMA_2P_GMR_filtered_zscored.Rds")
pyspi_data <- readRDS(pyspi_data_file) %>%
  mutate(group = stringr::str_to_sentence(group))

SPI_directionality <- read.csv(paste0(github_dir, "pairwise_analysis/SPI_Direction_Info.csv"))

################################################################################
# Generate model-free shuffle null distribution
################################################################################
if (!file.exists(paste0(rdata_path, sprintf("Null_Model_Free_Shuffles_%s.Rds",
                                            feature_set)))) {
  model_free_shuffle_null_res <- run_model_free_n_shuffles(num_shuffles = 100000,
                                                           feature_set = feature_set,
                                                           rdata_path = rdata_path)
  saveRDS(model_free_shuffle_null_res, file = paste0(rdata_path, sprintf("Null_Model_Free_Shuffles_%s.Rds",
                                                                         feature_set)))
} else {
  model_free_shuffle_null_res <- readRDS(paste0(rdata_path, sprintf("Null_Model_Free_Shuffles_%s.Rds",
                                                                    feature_set)))
}

################################################################################
# Create ten folds to use for all analyses
################################################################################

subjects_to_use <- readRDS(paste0(rdata_path, "UCLA_Subjects_with_Univariate_and_Pairwise.Rds"))

if (!file.exists(paste0(rdata_path, "Subjects_per_10_folds.Rds"))) {
  # Make folds
  set.seed(127)
  k = 10
  subject_folds <- caret::createFolds(subjects_to_use$group, k = k, list = TRUE, returnTrain = FALSE)
  
  # Save to Rds file
  saveRDS(subject_folds, file=paste0(rdata_path, "Subjects_per_10_folds.Rds"))
} else {
  subject_folds <- readRDS(paste0(rdata_path, "Subjects_per_10_folds.Rds"))
}

################################################################################
# Define weighting parameters
################################################################################

weighting_param_df <- data.frame(name = c("inv_prob"),
                                 use_inv_prob_weighting = c(TRUE),
                                 use_SMOTE = c(FALSE))

SVM_grouping_params <- data.frame(grouping_var = c("SPI"),
                                  SVM_feature_var = c("region_pair"))

################################################################################
# Run linear SVM for each grouping var
################################################################################

for (i in 1:nrow(weighting_param_df)) {
  weighting_name <- weighting_param_df$name[i]
  use_inv_prob_weighting <- weighting_param_df$use_inv_prob_weighting[i]
  use_SMOTE <- weighting_param_df$use_SMOTE[i]
  
  for (j in 1:nrow(SVM_grouping_params)) {
    grouping_var = SVM_grouping_params$grouping_var[j]
    SVM_feature_var = SVM_grouping_params$SVM_feature_var[j]
    
    # Run given weighting for 10-fold CV linear SVM
    if (!file.exists(paste0(rdata_path, sprintf("pyspi_%s_pairwise_CV_linear_SVM_%s_%s.Rds",
                                                grouping_var, feature_set, weighting_name)))) {
      tryCatch({group_wise_SVM_CV_weighting <- run_pairwise_cv_svm_by_input_var(pairwise_data = pyspi_data,
                                                                                           SPI_directionality = SPI_directionality,
                                                                                           svm_kernel = "linear",
                                                                                           num_k_folds = 10,
                                                                                           flds = subject_folds,
                                                                                           grouping_var = grouping_var,
                                                                                           svm_feature_var = SVM_feature_var,
                                                                                           test_package = "e1071",
                                                                                           noise_proc = "AROMA+2P+GMR",
                                                                                           use_inv_prob_weighting = use_inv_prob_weighting,
                                                                                           use_SMOTE = use_SMOTE,
                                                                                           shuffle_labels = FALSE)
      saveRDS(group_wise_SVM_CV_weighting, file=paste0(rdata_path,
                                                                  sprintf("pyspi_%s_pairwise_CV_linear_SVM_%s_%s.Rds",
                                                                          grouping_var,
                                                                          feature_set,
                                                                          weighting_name)))
      }, error = function(e) {
        cat("\nCould not run", grouping_var, "wise analysis:\n")
        print(e)
      })
    } else {
      group_wise_SVM_CV_weighting <- readRDS(paste0(rdata_path,
                                                    sprintf("pyspi_%s_pairwise_CV_linear_SVM_%s_%s.Rds",
                                                            grouping_var,
                                                            feature_set,
                                                            weighting_name)))
    }
    
    
    #### Calculate balanced accuracy across all folds
    if (!file.exists(paste0(rdata_path, sprintf("pyspi_%s_pairwise_CV_linear_SVM_%s_%s_balacc.Rds",
                                                grouping_var, feature_set, weighting_name)))) {
      group_wise_SVM_balanced_accuracy <- group_wise_SVM_CV_weighting %>%
        group_by(grouping_var, Noise_Proc, Sample_Type) %>%
        summarise(accuracy = sum(Prediction_Correct) / n(),
                  balanced_accuracy = caret::confusionMatrix(data = Predicted_Diagnosis,
                                                             reference = Actual_Diagnosis)$byClass[["Balanced Accuracy"]])
      
      saveRDS(group_wise_SVM_balanced_accuracy, file=paste0(rdata_path, sprintf("pyspi_%s_pairwise_CV_linear_SVM_%s_%s_balacc.Rds",
                                                                                grouping_var, feature_set, weighting_name)))
    }
    
    #### Calculate p values from model-free shuffle null distribution
    if (!file.exists(paste0(rdata_path, sprintf("pyspi_%s_pairwise_CV_linear_SVM_model_free_shuffle_pvals_%s_%s.Rds",
                                                grouping_var, feature_set, weighting_name)))) {
      group_wise_SVM_balanced_accuracy <- readRDS(paste0(rdata_path,
                                                        sprintf("pyspi_%s_pairwise_CV_linear_SVM_%s_%s_balacc.Rds",
                                                                grouping_var,
                                                                feature_set, 
                                                                weighting_name)))
      
      # Calculate p-values
      pvalues <- calc_empirical_nulls(class_res = group_wise_SVM_balanced_accuracy,
                                      null_data = model_free_shuffle_null_res,
                                      feature_set = feature_set,
                                      use_pooled_null = TRUE,
                                      is_main_data_averaged = TRUE,
                                      grouping_var = grouping_var)
      
      saveRDS(pvalues, file=paste0(rdata_path, sprintf("pyspi_%s_pairwise_CV_linear_SVM_model_free_shuffle_pvals_%s_%s.Rds",
                                                       grouping_var, feature_set, weighting_name)))
    }
    
    weighting_null_dist_file <- paste0(rdata_path, sprintf("pyspi_%s_pairwise_%s_%s_null_model_fits.Rds",
                                                           grouping_var, feature_set, weighting_name))
    
    # Run null perm iterations if overall null distribution data file doesn't exist
    if (!file.exists(weighting_null_dist_file)) {
      # Output script dir
      output_data_dir <- paste0(rdata_path, sprintf("pyspi_%s_pairwise_%s_%s_null_model_fits/",
                                                    grouping_var, feature_set, weighting_name))
      output_scripts_dir <- paste0(github_dir, sprintf("pairwise_analysis/pyspi_%s_pairwise_%s_%s_null_model_fits/",
                                                       grouping_var, weighting_name, feature_set))
      icesTAF::mkdir(output_scripts_dir)
      # template file for null distributions
      if (grouping_var == "SPI") {
        num_permutations <- 100
        nperm_per_iter <- 10
      } else {
        num_permutations <- 10
        nperm_per_iter <- 1
      }
      num_k_folds <- 10
      template_pbs_file <- paste0(github_dir, "pairwise_analysis/template_null_model_fit.pbs")
      
      
      lookup_list <- list("PROJECT_NAME" = "hctsa",
                          "NAME" = sprintf("pyspi_%swise_null_model_fit", grouping_var),
                          "MEMNUM" = "20",
                          "NCPUS" = "1",
                          "GITHUB_DIR" = github_dir,
                          "PROJECT_DIR" = project_path,
                          "EMAIL" = "abry4213@uni.sydney.edu.au",
                          "PBS_NOTIFY" = "a",
                          "WALL_HRS" = "6",
                          "PAIRWISE_DATA_FILE" = paste0(pydata_path, sprintf("UCLA_all_subject_%s_AROMA_2P_GMR_filtered_zscored.Rds",
                                                                             feature_set)),
                          "SPI_DIRECTIONALITY_FILE" = paste0(github_dir, "pairwise_analysis/SPI_Direction_Info.csv"),
                          "NUM_K_FOLDS" = num_k_folds,
                          "NUM_PERMS_PER_ITER" = nperm_per_iter,
                          "OUTPUT_DATA_DIR" = output_data_dir,
                          "FEATURE_SET" = "pyspi_19",
                          "GROUPING_VAR" = grouping_var,
                          "SVM_FEATURE_VAR" = SVM_feature_var,
                          "NOISE_PROC" = noise_proc)
      
      to_be_replaced <- names(lookup_list)
      replacement_values <- unlist(unname(lookup_list))
      
      for (p in 1:num_permutations) {
        
        
        # Run command if null file doesn't exist
        if (!file.exists(sprintf("%s/pyspi_%s_pairwise_%s_inv_prob_null_model_fit_iter_%s.Rds",
                                 output_data_dir, grouping_var, feature_set, p))) {
          cat("\nNow running null perms for iteration", p, "\n")
          new_pbs_file <- readLines(template_pbs_file)
          
          # Replace file paths
          pbs_text_replaced <- mgsub::mgsub(new_pbs_file,
                                            to_be_replaced,
                                            replacement_values)
          
          # Replace null iteration number
          pbs_text_replaced <- gsub("iterj", p, pbs_text_replaced)
          
          # Write updated PBS script to file
          output_pbs_file <- writeLines(pbs_text_replaced,
                                        paste0(output_scripts_dir,
                                               "null_iter_", p, ".pbs"))
          
          system(paste0("qsub ", output_scripts_dir, "null_iter_", p, ".pbs"))
          
        }
      }
    }
  }
}