################################################################################
# Load libraries
################################################################################

library(tidyverse)
library(icesTAF)
library(cowplot)
library(ggpubr)
library(ggsignif)
theme_set(theme_cowplot())

################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
source(paste0(github_dir, "helper_functions/classification/Linear_SVM.R"))
source(paste0(github_dir, "data_visualisation/manuscript_figures/Manuscript_Draft_Visualisations_Helper.R"))
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/FigureS1/")
TAF::mkdir(plot_path)

SCZ_data_path <- "~/data/UCLA_Schizophrenia/"
SCZ_rdata_path <- paste0(SCZ_data_path, "processed_data/Rdata/")

# Load subject metadata
SCZ_subject_metadata <- readRDS(paste0(SCZ_data_path, "UCLA_Schizophrenia_sample_metadata.Rds"))

# Load fractional displacement (movement) data
SCZ_movement_data <- compile_movement_data(fd_path = paste0(SCZ_data_path, "movementData/"),
                                           input_dataset_name = "UCLA_Schizophrenia",
                                           sample_metadata = SCZ_subject_metadata) %>%
  filter(!is.na(Diagnosis))

################################################################################
# Plot subject movement by diagnosis group
################################################################################

wilcox.test(FD ~ Diagnosis, data = SCZ_movement_data)

SCZ_movement_data %>%
  ggplot(data=., mapping=aes(x=Diagnosis, y=FD)) +
  geom_violin(aes(fill=Diagnosis)) +
  geom_boxplot(color="black", fill=NA, width=0.1) +
  ylab("Fractional Displacement (FD)") +
  xlab("Group") +
  scale_fill_manual(values = c("#00B06D", "#737373")) +
  theme(legend.position="none",
        plot.title=element_text(hjust=0.5)) +
  geom_signif(test = "wilcox.test",
              comparisons = list(c("Schizophrenia", "Control")), 
              map_signif_level=TRUE)
ggsave(paste0(plot_path, "SCZ_FD_by_group_violin.png"),
       width = 4, height = 4, units="in", dpi=300, bg="white")

################################################################################
# Plot # subjects retained per group by FD threshold
################################################################################
fd_thresh_list <- list()
# Iterate over thresholds from 0 to 1 at intervals of 0.05
for (fd_threshold in seq(0, 1, by=0.05)) {
  # Find number of subjects retained at this threshold per group
  num_samples_by_dx <- table(subset(SCZ_movement_data, FD <= fd_threshold)$Diagnosis)
  tryCatch({
    thresh_df <<- as.data.frame(num_samples_by_dx) %>%
      mutate(FD_Threshold = fd_threshold) %>%
      dplyr::rename("Diagnosis" = "Var1",
                    "Num_Subjects" = "Freq")
  }, error = function(e) {
    thresh_df <<- data.frame(Diagnosis = unique(SCZ_movement_data$Diagnosis),
                            Num_Subjects = 0,
                            FD_Threshold = fd_threshold)
  })

  fd_thresh_list <- rlist::list.append(fd_thresh_list, thresh_df)
}
SCZ_threshold_data <- do.call(plyr::rbind.fill, fd_thresh_list)

# Plot # of subjects retained per FD threshold by group
SCZ_threshold_data %>%
  group_by(Diagnosis) %>%
  ggplot(data=., mapping=aes(x=FD_Threshold, y=Num_Subjects, 
                             color=Diagnosis, group=Diagnosis)) +
  geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=0.12, xmax=0.5),
            color=NA, fill="gray85", alpha=0.5) +
  geom_line(size=1.75, alpha=0.9) +
  ylab("# Subjects") +
  xlab("Fractional Displacement (FD) Maximum Threshold") +
  scale_color_manual(values = c("#00B06D", "#737373")) +
  scale_x_reverse() +
  theme(legend.position="none",
        strip.placement = "outside",
        plot.title=element_text(hjust=0.5))
ggsave(paste0(plot_path, "SCZ_FD_threshold_num_subjects.png"),
       width = 5, height = 2.75, units="in", dpi=300, bg="white")

# Plot % of subjects retained per FD threshold by group
SCZ_threshold_data %>%
  group_by(Diagnosis) %>%
  mutate(Perc_Subjects = Num_Subjects / max(Num_Subjects, na.rm=T)) %>%
  ggplot(data=., mapping=aes(x=FD_Threshold, y=Perc_Subjects, 
                             color=Diagnosis, group=Diagnosis)) +
  geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=0.12, xmax=0.5),
            color=NA, fill="gray85", alpha=0.5) +
  geom_line(size=1.75, alpha=0.9) +
  ylab("% of Subjects") +
  xlab("Fractional Displacement (FD) Maximum Threshold") +
  scale_color_manual(values = c("#00B06D", "#737373")) +
  scale_x_reverse() +
  theme(legend.position="bottom",
        strip.placement = "outside",
        plot.title=element_text(hjust=0.5))
ggsave(paste0(plot_path, "SCZ_FD_threshold_percent_subjects.png"),
       width = 5, height = 3, units="in", dpi=300, bg="white")

################################################################################
# Balanced accuracy as a function of FD threshold for univariate combo
################################################################################

# Define constants
grouping_type <- "Combo"
grouping_var <- "Combo"
SVM_feature_var <- "Combo"
top_region <- "ctx-rh-postcentral"
univariate_feature_set <- "catch22"
pairwise_feature_set <- "pyspi14"
noise_proc <- "AROMA+2P+GMR"
num_k_folds <- 10
nrepeats <- 10
svm_kernel <- "linear"

# Load SCZ catch22 z-scored data
SCZ_catch22_zscored <- readRDS(paste0(SCZ_rdata_path, 
                                      "UCLA_Schizophrenia_catch22_filtered_zscored.Rds")) %>%
  filter(Noise_Proc == noise_proc)


# Get diagnosis proportions
SCZ_sample_groups <- readRDS(paste0(SCZ_rdata_path, 
                                sprintf("UCLA_Schizophrenia_samples_with_univariate_%s_and_pairwise_%s_filtered.Rds",
                                                    univariate_feature_set,
                                                    pairwise_feature_set))) %>%
  left_join(., SCZ_subject_metadata) %>%
  distinct(Sample_ID, Diagnosis)

# Function to run 10-repeat 10-fold linear SVM 
run_repeat_cv_linear_svm <- function(fd_list = seq(0.12, 0.5, by=0.02),
                                     movement_data,
                                     catch22_data,
                                     type = "Brain Region",
                                     input_region = "") {
  # Iterate over each threshold
  fd_thresh_list <- list()
  
  # Iterate over thresholds from 0 to 1 at intervals of 0.05
  for (fd_threshold in seq(0.12, 0.5, by=0.02)) {
    # Data thresholded by FD
    movement_data_thresh <- subset(movement_data, FD <= fd_threshold)
    
    catch22_data_thresh <- subset(catch22_data, 
                                  Sample_ID %in% movement_data_thresh$Sample_ID)
    
    # Define sample weights for inverse probability weighting
    sample_wts <- as.list(1/prop.table(table(movement_data_thresh$Diagnosis)))
    
    if (type=="Brain Region") {
      # Prep data for SVM
      data_for_SVM <- catch22_data_thresh %>%
        filter(Brain_Region == input_region) %>%
        dplyr::ungroup() %>%
        left_join(., SCZ_sample_groups) %>%
        dplyr::select(Sample_ID, Diagnosis, Brain_Region, values) %>%
        distinct(.keep_all = T) %>%
        tidyr::pivot_wider(id_cols = c(Sample_ID, Diagnosis),
                           names_from = Brain_Region,
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
        left_join(., SCZ_sample_groups) %>%
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
    } else if (type=="FD Only") {
      type_label = "FD Only"
      data_for_SVM <- movement_data_thresh %>%
        dplyr::select(Sample_ID, Diagnosis, FD)
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
                                       FD_threshold = fd_threshold))
      fd_thresh_list <- rlist::list.append(fd_thresh_list, SVM_results)
    }
  }
  SVM_res <- do.call(plyr::rbind.fill, fd_thresh_list)  %>%
    group_by(FD_threshold, repeat_number) %>%
    summarise(balanced_accuracy = caret::confusionMatrix(data = Predicted_Diagnosis,
                                                         reference = Actual_Diagnosis)$byClass[["Balanced Accuracy"]]) %>%
    group_by(FD_threshold) %>%
    summarise(meanbacc = 100*mean(balanced_accuracy),
              sdbacc = 100*sd(balanced_accuracy)) %>%
    mutate(Method = type_label)
}

combo_catch22_svm_res <- run_repeat_cv_linear_svm(movement_data = SCZ_movement_data,
                         catch22_data = SCZ_catch22_zscored,
                         type = "Combo"
)
saveRDS(SCZ_threshold_SVM_res, "SCZ_threshold_SVM_res.Rds")

SCZ_threshold_SVM_res %>%
  plyr::rbind.fill(., SCZ_threshold_SVM_mvmt_res) %>%
  ggplot(data=., mapping=aes(x=FD_threshold)) +
  scale_x_reverse() +
  geom_ribbon(aes(ymin = meanbacc - sdbacc,
                  ymax = meanbacc + sdbacc,
                  fill = Method),
              alpha=0.2) +
  geom_line(aes(y=meanbacc, color = Method)) +
  xlab("FD Maximum Threshold") +
  ylab("Balanced Accuracy (%)") 
ggsave(paste0(plot_path, "SCZ_FD_threshold_SVM_balanced_accuracy.png"),
       width = 5, height = 2.5, units="in", dpi=300, bg="white")
