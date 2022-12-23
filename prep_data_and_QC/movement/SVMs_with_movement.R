################################################################################
# Load libraries
################################################################################

library(tidyverse)
library(icesTAF)

################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
source(paste0(github_dir, "helper_functions/classification/Linear_SVM.R"))
source(paste0(github_dir, "data_visualisation/manuscript_figures/Manuscript_Draft_Visualisations_Helper.R"))
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/FigureS1/")
TAF::mkdir(plot_path)

data_path <- "~/data/UCLA_CNP_ABIDE_ASD/"
rdata_path <- paste0(data_path, "processed_data/Rdata/")

univariate_feature_set <- "catch22"
pairwise_feature_set <- "pyspi14"
num_k_folds <- 10
nrepeats <- 10
svm_kernel <- "linear"

################################################################################
# Load data
################################################################################

# Constants
UCLA_top_region <- "ctx-rh-postcentral"
UCLA_noise_proc <- "AROMA+2P+GMR"
ABIDE_noise_proc <- "FC1000"
ABIDE_top_region <- "Superior Frontal Gyrus"

# Load subject metadata
sample_metadata <- readRDS(paste0(data_path, "study_metadata/UCLA_CNP_ABIDE_ASD_sample_metadata.Rds"))

if (!file.exists(paste0(data_path, "movement_data/UCLA_CNP_ABIDE_ASD_FD.Rds"))) {
  # Load Annie's mean FD (Jenkinson, Power, VanDijk) data
  UCLA_movement_data <- read.table(paste0(data_path, "movement_data/UCLA_CNP/UCLA_CNP_mFD.txt"),
                                   sep=",")
  colnames(UCLA_movement_data) <- c("Sample_ID", "Jenkinson", "Power", "VanDijk")
  UCLA_movement_data <- left_join(UCLA_movement_data, sample_metadata) 
  
  
  # Load Linden's mean FD (Power) data
  UCLA_movement_data_Linden <- read.table(paste0(data_path, 
                                                 "movement_data/UCLA_CNP/fdAvgs_UCLA_Schizophrenia.txt"),
                                          sep=",")
  colnames(UCLA_movement_data_Linden) <- c("Jenkinson_Linden")
  UCLA_movement_data_Linden <- cbind(UCLA_movement_data_Linden, UCLA_movement_data)
  # Save the data with Linden's estimates
  saveRDS(UCLA_movement_data_Linden, file=paste0(data_path, "movement_data/UCLA_CNP_Linden_FD_Power.Rds"))
  
  # Filter data down for UCLA
  UCLA_movement_data <- UCLA_movement_data %>%
    filter(!is.na(Diagnosis)) %>%
    mutate(Cohort="UCLA Study")
  
  # Load mean framewise displacement (movement) data for ABIDE ASD + Control participants
  ABIDE_movement_data <- read.table(paste0(data_path, "movement_data/ABIDE_ASD/ABIDE_ASD_mFD.txt"),
                                    sep=",", colClasses = c("V1" = "character"))
  colnames(ABIDE_movement_data) <- c("Sample_ID", "Jenkinson", "Power", "VanDijk")
  ABIDE_movement_data <- left_join(ABIDE_movement_data, sample_metadata) %>%
    filter(!is.na(Diagnosis)) %>%
    mutate(Cohort="ABIDE Study")
  
  # Combine data and save
  UCLA_ABIDE_movement_data <- plyr::rbind.fill(UCLA_movement_data, ABIDE_movement_data)
  saveRDS(UCLA_ABIDE_movement_data, file=paste0(data_path, "movement_data/UCLA_CNP_ABIDE_ASD_FD.Rds"))
} else {
  UCLA_ABIDE_movement_data <- readRDS(paste0(data_path, "movement_data/UCLA_CNP_ABIDE_ASD_FD.Rds"))
}

# Load catch22 z-scored data
UCLA_catch22_zscored <- readRDS(paste0(rdata_path, 
                                       "UCLA_CNP_catch22_filtered_zscored.Rds")) %>%
  left_join(., sample_metadata) %>%
  filter(Noise_Proc == UCLA_noise_proc,
         Diagnosis %in% c("Control", "Schizophrenia"))

ABIDE_catch22_zscored <- readRDS(paste0(rdata_path, 
                                        "ABIDE_ASD_catch22_filtered_zscored.Rds")) %>%
  filter(Noise_Proc == ABIDE_noise_proc)

# Get diagnosis proportions
UCLA_sample_groups <- data.frame(Sample_ID = readRDS(paste0(rdata_path, 
                                                            sprintf("UCLA_CNP_ABIDE_ASD_samples_with_univariate_%s_and_pairwise_%s_filtered.Rds",
                                                                    univariate_feature_set,
                                                                    pairwise_feature_set)))) %>%
  left_join(., sample_metadata) %>%
  filter(Diagnosis %in% c("Control", "Schizophrenia"),
         Study == "UCLA_CNP") %>%
  distinct(Sample_ID, Diagnosis)

ABIDE_sample_groups <- readRDS(paste0(rdata_path, 
                                      sprintf("ABIDE_ASD_samples_with_univariate_%s_and_pairwise_%s_filtered.Rds",
                                              univariate_feature_set,
                                              pairwise_feature_set))) %>%
  left_join(., sample_metadata) %>%
  distinct(Sample_ID, Diagnosis)

# Define movement ranges per study
UCLA_seq_list <- seq(0.12, max(subset(UCLA_ABIDE_movement_data,
                                      Study=="UCLA_CNP") %>%
                                 pull(Power)), by=0.02)
ABIDE_seq_list <- seq(0.01, 1, by=0.04)

################################################################################
# Balanced accuracy as a function of FD threshold for univariate combo
################################################################################

# Function to run 10-repeat 10-fold linear SVM 
run_repeat_cv_linear_svm <- function(mvmt_list = seq(0.12, 1, by=0.05),
                                     movement_data,
                                     movement_var = "Power",
                                     sample_groups,
                                     catch22_data,
                                     type = "Brain Region",
                                     input_region = "") {
  # Iterate over each threshold
  df_list <- list()
  
  # Iterate over thresholds from 0 to 1 at intervals of 0.05
  for (movement_threshold in mvmt_list) {
    # Data thresholded by FD
    movement_data_thresh <- movement_data %>%
      mutate(movement_var = get(movement_var)) %>%
      filter(movement_var <= movement_threshold)
    
    catch22_data_thresh <- subset(catch22_data, 
                                  Sample_ID %in% movement_data_thresh$Sample_ID)
    
    # Define sample weights for inverse probability weighting
    sample_wts <- as.list(1/prop.table(table(movement_data_thresh$Diagnosis)))
    
    if (type=="Brain Region") {
      # Prep data for SVM
      data_for_SVM <- catch22_data_thresh %>%
        filter(Brain_Region == input_region) %>%
        dplyr::ungroup() %>%
        left_join(., sample_groups) %>%
        dplyr::select(Sample_ID, Diagnosis, names, values) %>%
        distinct(.keep_all = T) %>%
        tidyr::pivot_wider(id_cols = c(Sample_ID, Diagnosis),
                           names_from = names,
                           values_from 
                           = values) %>%
        # Drop columns that are all NA/NAN
        dplyr::select(where(function(x) any(!is.na(x)))) %>%
        # Drop rows with NA for one or more column
        drop_na()
      
      type_label = input_region
      
    } else if (type=="Combo") {
      type_label = "catch22 Combo"
      # Prep data for SVM
      data_for_SVM <- catch22_data_thresh %>%
        unite("Combo", c("Brain_Region", "names"), sep="_", remove=F) %>%
        dplyr::ungroup() %>%
        left_join(., sample_groups) %>%
        dplyr::select(Sample_ID, Diagnosis, Combo, values) %>%
        distinct(.keep_all = T) %>%
        tidyr::pivot_wider(id_cols = c(Sample_ID, Diagnosis),
                           names_from = Combo,
                           values_from 
                           = values) %>%
        # Drop columns that are all NA/NAN
        dplyr::select(where(function(x) any(!is.na(x)))) %>%
        # Drop rows with NA for one or more column
        drop_na()
    } else if (type=="Movement Only") {
      type_label = "Movement Only"
      data_for_SVM <- movement_data_thresh %>%
        dplyr::select(Sample_ID, Diagnosis, movement_var)
    }
    
    # Run linear SVM
    if (nrow(data_for_SVM) > 0) {
      SVM_results <- 1:nrepeats %>%
        purrr::map_df( ~ k_fold_CV_linear_SVM(input_data = data_for_SVM,
                                              k = num_k_folds,
                                              svm_kernel = svm_kernel,
                                              sample_wts = sample_wts,
                                              shuffle_labels = F,
                                              out_of_sample_only = T)%>%
                         dplyr::mutate(repeat_number = .x,
                                       movement_threshold = movement_threshold))
      df_list <- list.append(df_list, SVM_results)
    }
  }
  SVM_res <- do.call(plyr::rbind.fill, df_list)  %>%
    group_by(movement_threshold, repeat_number) %>%
    summarise(balanced_accuracy = caret::confusionMatrix(data = Predicted_Diagnosis,
                                                         reference = Actual_Diagnosis)$byClass[["Balanced Accuracy"]]) %>%
    group_by(movement_threshold) %>%
    summarise(meanbacc = 100*mean(balanced_accuracy),
              sdbacc = 100*sd(balanced_accuracy)) %>%
    mutate(Method = type_label)
  
  return(SVM_res)
}

# Run for UCLA CNP
if (!file.exists(paste0(rdata_path, "UCLA_CNP_catch22_SVM_with_Movement_Thresholds.Rdata"))) {

    # Run SVM with various FD threshold cutoffs using just movement data
  UCLA_mvmt_svm_res <- run_repeat_cv_linear_svm(movement_data = subset(UCLA_ABIDE_movement_data, Study=="UCLA_CNP"),
                                                mvmt_list = UCLA_seq_list,
                                                catch22_data = UCLA_catch22_zscored,
                                                movement_var = "Power",
                                                type = "Movement Only",
                                                sample_groups = UCLA_sample_groups
  )

  # Run SVM with various FD threshold cutoffs using right postcentral cortex catch22
  UCLA_ROI_catch22_svm_res <- run_repeat_cv_linear_svm(movement_data = subset(UCLA_ABIDE_movement_data, Study=="UCLA_CNP"),
                                                      mvmt_list = UCLA_seq_list,
                                                      catch22_data = UCLA_catch22_zscored,
                                                      movement_var = "Power",
                                                      type = "Brain Region",
                                                      input_region = UCLA_top_region,
                                                      sample_groups = UCLA_sample_groups
  )

  # Run SVM with various FD threshold cutoffs using univariate combo catch22
  UCLA_combo_catch22_svm_res <- run_repeat_cv_linear_svm(
    mvmt_list = UCLA_seq_list,
    movement_data = subset(UCLA_ABIDE_movement_data, Study=="UCLA_CNP"),
    catch22_data = UCLA_catch22_zscored,
    movement_var = "Power",
    type = "Combo",
    sample_groups = UCLA_sample_groups
  )
  save(UCLA_combo_catch22_svm_res,
      UCLA_ROI_catch22_svm_res,
      UCLA_mvmt_svm_res,
      file=paste0(rdata_path, "UCLA_CNP_catch22_SVM_with_Movement_Thresholds.Rdata"))

}
# Run for ABIDE ASD
if (!file.exists(paste0(rdata_path, "ABIDE_ASD_catch22_SVM_with_Movement_Thresholds.Rdata"))) {
  
  # Run SVM with various FD threshold cutoffs using just movement data
  ABIDE_mvmt_svm_res <- run_repeat_cv_linear_svm(
    mvmt_list = ABIDE_seq_list,
    movement_data = subset(UCLA_ABIDE_movement_data, Study=="ABIDE_ASD"),
    catch22_data = ABIDE_catch22_zscored,
    movement_var = "Power",
    type = "Movement Only",
    sample_groups = ABIDE_sample_groups
  )
  # Run SVM with various FD threshold cutoffs using superior frontal gyrus catch22
  ABIDE_ROI_catch22_svm_res <- run_repeat_cv_linear_svm(
    mvmt_list = ABIDE_seq_list,
    movement_data = subset(UCLA_ABIDE_movement_data, Study=="ABIDE_ASD"),
    catch22_data = ABIDE_catch22_zscored,
    movement_var = "Power",
    type = "Brain Region",
    input_region = ABIDE_top_region,
    sample_groups = ABIDE_sample_groups
  )
  # Run SVM with various FD threshold cutoffs using univariate combo catch22
  ABIDE_combo_catch22_svm_res <- run_repeat_cv_linear_svm(
    mvmt_list = ABIDE_seq_list,
    movement_data = subset(UCLA_ABIDE_movement_data, Study=="ABIDE_ASD"),
    catch22_data = ABIDE_catch22_zscored,
    movement_var = "Power",
    type = "Combo",
    sample_groups = ABIDE_sample_groups
  )
  save(ABIDE_combo_catch22_svm_res,
     ABIDE_ROI_catch22_svm_res,
     ABIDE_mvmt_svm_res,
     file=paste0(rdata_path, "ABIDE_ASD_catch22_SVM_with_Movement_Thresholds.Rdata"))
}



