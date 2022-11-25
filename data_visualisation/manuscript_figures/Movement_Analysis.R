################################################################################
# Load libraries
################################################################################

library(tidyverse)
library(icesTAF)
library(cowplot)
library(ggpubr)
library(ggsignif)
library(patchwork)
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
ASD_data_path <- "~/data/ABIDE_ASD/"
ASD_rdata_path <- paste0(ASD_data_path, "processed_data/Rdata/")

# Load subject metadata
SCZ_subject_metadata <- readRDS(paste0(SCZ_data_path, "UCLA_Schizophrenia_sample_metadata.Rds"))
ASD_subject_metadata <- readRDS(paste0(ASD_data_path, "ABIDE_ASD_sample_metadata.Rds"))

# Load fractional displacement (movement) data for UCLA Schizophrenia
SCZ_movement_data <- compile_movement_data(fd_path = paste0(SCZ_data_path, "movementData/"),
                                           input_dataset_name = "UCLA_Schizophrenia",
                                           sample_metadata = SCZ_subject_metadata) %>%
  filter(!is.na(Diagnosis))

# Load mean displacement data for ABIDE ASD
ASD_motion_path <- paste0(ASD_data_path, "raw_data/movement_data/")
ASD_movement_data <- calculate_mean_displacement(movement_data_path = ASD_motion_path,
                                                   input_dataset_name = "ABIDE_ASD",
                                                   sample_metadata = ASD_subject_metadata)

################################################################################
# Plot subject movement by diagnosis group
################################################################################

# ABIDE ASD
wilcox.test(mean_displacement ~ Diagnosis, data = ASD_movement_data)

ASD_p <- ASD_movement_data %>%
  filter(mean_displacement<=5) %>%
  mutate(Diagnosis = factor(Diagnosis, levels = c("Control", "ASD")),
         Method = "ASD") %>%
  ggplot(data=., mapping=aes(x=Diagnosis, y=mean_displacement)) +
  geom_violin(aes(fill=Diagnosis)) +
  geom_boxplot(color="black", fill=NA, width=0.1) +
  scale_fill_manual(values = c("#00B06D", "#737373")) +
  scale_y_continuous(expand = c(0,0,0.15,0)) +
  xlab("Group") +
  ylab("Mean\nDisplacement (mm)") +
  facet_grid(. ~ Method) +
  theme(legend.position="none",
        plot.title=element_text(hjust=0.5)) +
  geom_signif(test = "wilcox.test",
              comparisons = list(c("ASD", "Control")), 
              map_signif_level=TRUE)

# UCLA Schizophrenia
wilcox.test(FD ~ Diagnosis, data = SCZ_movement_data)

SCZ_p <- SCZ_movement_data %>%
  mutate(Diagnosis = factor(Diagnosis, levels = c("Control", "Schizophrenia")),
         Method = "SCZ") %>%
  ggplot(data=., mapping=aes(x=Diagnosis, y=FD)) +
  geom_violin(aes(fill=Diagnosis)) +
  geom_boxplot(color="black", fill=NA, width=0.1) +
  xlab("Group") +
  ylab("Fractional\nDisplacement") +
  scale_fill_manual(values = c("#00B06D", "#737373")) +
  scale_y_continuous(expand = c(0,0,0.15,0)) +
  facet_grid(. ~ Method) +
  theme(legend.position="none",
        plot.title=element_text(hjust=0.5)) +
  geom_signif(test = "wilcox.test",
              comparisons = list(c("Schizophrenia", "Control")), 
              map_signif_level=TRUE)

ASD_p+SCZ_p
ggsave(paste0(plot_path, "ASD_SCZ_Movement_by_Group.png"),
       width = 6, height=2.5, units="in", dpi=300)

################################################################################
# Plot # subjects retained per group by motion threshold
################################################################################
# Function to find number of subjects retained per group by motion threshold
subjects_retained_by_motion <- function(movement_data, 
                                        motion_variable = "FD",
                                        motion_range=seq(0, 1, by=0.05)) {
  df_list <- list()
  # Iterate over thresholds from 0 to 1 at intervals of 0.05
  for (motion_threshold in motion_range) {
    # Find number of subjects retained at this threshold per group
    num_samples_by_dx <- table(movement_data %>%
                                 filter(get(motion_variable) <= motion_threshold) %>%
                                 pull(Diagnosis))
    tryCatch({
      mvmt_data <<- as.data.frame(num_samples_by_dx) %>%
        mutate(Motion_Threshold = motion_threshold) %>%
        dplyr::rename("Diagnosis" = "Var1",
                      "Num_Subjects" = "Freq")
    }, error = function(e) {
      mvmt_data <<- data.frame(Diagnosis = unique(movement_data$Diagnosis),
                               Num_Subjects = 0,
                               Motion_Threshold = motion_threshold)
    })
    
    df_list <- list.append(df_list, mvmt_data)
  }
  results_df <- do.call(plyr::rbind.fill, df_list)
  return(results_df)
}

# Find # subjects retained per motion threshold by group
SCZ_subjects_by_motion <- subjects_retained_by_motion(movement_data = SCZ_movement_data,
                                                      motion_variable = "FD",
                                                      motion_range=seq(0, 1, by=0.05))
ASD_subjects_by_motion <- subjects_retained_by_motion(movement_data = ASD_movement_data,
                                                      motion_variable = "mean_displacement",
                                                      motion_range=c(seq(0,0.5, by=0.005),
                                                                     seq(0.5, 5, by=0.25)))


# Plot # of subjects retained per FD threshold by group
# SCZ
SCZ_num_motion <- SCZ_subjects_by_motion %>%
  group_by(Diagnosis) %>%
  ggplot(data=., mapping=aes(x=Motion_Threshold, y=Num_Subjects, 
                             color=Diagnosis, group=Diagnosis)) +
  geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=0.12, xmax=0.5),
            color=NA, fill="gray85", alpha=0.5) +
  geom_line(size=1.75, alpha=0.9) +
  ylab("# Subjects") +
  scale_color_manual(values = c("#00B06D", "#737373")) +
  scale_x_reverse() +
  theme(legend.position="none",
        axis.title.x = element_blank(),
        strip.placement = "outside",
        plot.title=element_text(hjust=0.5))
# ASD
ASD_num_motion <- ASD_subjects_by_motion %>%
  mutate(Diagnosis = factor(Diagnosis, levels=c("Control", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Motion_Threshold, y=Num_Subjects, 
                             color=Diagnosis, group=Diagnosis)) +
  geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=0.05, xmax=0.6),
            color=NA, fill="gray85", alpha=0.5) +
  geom_line(size=1.75, alpha=0.8) +
  ylab("# Subjects") +
  scale_color_manual(values = c("#00B06D", "#737373")) +
  scale_x_reverse() +
  theme(legend.position="none",
        axis.title.x = element_blank(),
        strip.placement = "outside",
        plot.title=element_text(hjust=0.5))

# Plot % of subjects retained per FD threshold by group
# SCZ
SCZ_perc_motion <- SCZ_subjects_by_motion %>%
  group_by(Diagnosis) %>%
  mutate(Perc_Subjects = Num_Subjects / max(Num_Subjects, na.rm=T)) %>%
  ggplot(data=., mapping=aes(x=Motion_Threshold, y=Perc_Subjects, 
                             color=Diagnosis, group=Diagnosis)) +
  geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=0.12, xmax=0.5),
            color=NA, fill="gray85", alpha=0.5) +
  geom_line(size=1.75, alpha=0.9) +
  ylab("% of Subjects") +
  labs(color="Group") +
  xlab("Fractional Displacement (FD)\nMaximum Threshold") +
  scale_color_manual(values = c("#00B06D", "#737373")) +
  scale_x_reverse() +
  theme(legend.position="bottom",
        strip.placement = "outside",
        plot.title=element_text(hjust=0.5))
# ASD
ASD_perc_motion <- ASD_subjects_by_motion %>%
  mutate(Diagnosis = factor(Diagnosis, levels=c("Control", "ASD"))) %>%
  group_by(Diagnosis) %>%
  mutate(Perc_Subjects = Num_Subjects / max(Num_Subjects, na.rm=T)) %>%
  ggplot(data=., mapping=aes(x=Motion_Threshold, y=Perc_Subjects, 
                             color=Diagnosis, group=Diagnosis)) +
  geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=0.05, xmax=0.6),
            color=NA, fill="gray85", alpha=0.5) +
  geom_line(size=1.75, alpha=0.9) +
  ylab("% of Subjects") +
  labs(color="Group") +
  xlab("Mean Displacement (mm)\nMaximum Threshold") +
  scale_color_manual(values = c("#00B06D", "#737373")) +
  scale_x_reverse() +
  theme(legend.position="bottom",
        strip.placement = "outside",
        plot.title=element_text(hjust=0.5))

ASD_num_motion/ASD_perc_motion
ggsave(paste0(plot_path, "ASD_mvmt_threshold_num_subjects.png"),
       width = 5, height = 4, units="in", dpi=300, bg="white")

SCZ_num_motion/SCZ_perc_motion
ggsave(paste0(plot_path, "SCZ_mvmt_threshold_num_subjects.png"),
       width = 5, height = 4, units="in", dpi=300, bg="white")

################################################################################
# Balanced accuracy as a function of FD threshold for univariate combo
################################################################################

univariate_feature_set <- "catch22"
pairwise_feature_set <- "pyspi14"
num_k_folds <- 10
nrepeats <- 10
svm_kernel <- "linear"

# Function to run 10-repeat 10-fold linear SVM 
run_repeat_cv_linear_svm <- function(mvmt_list = seq(0.12, 0.5, by=0.02),
                                     movement_data,
                                     movement_var = "FD",
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
      tryCatch({
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
      }, error = function(e) {
        
      })
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
}

# Define constants for UCLA Schizophrenia
top_region <- "ctx-rh-postcentral"
noise_proc <- "AROMA+2P+GMR"

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


# Run SVM with various FD threshold cutoffs using univariate combo catch22
SCZ_combo_catch22_svm_res <- run_repeat_cv_linear_svm(
  mvmt_list = seq(0.12, 0.5, by=0.02),
  movement_data = SCZ_movement_data,
  catch22_data = SCZ_catch22_zscored,
  movement_var = "FD",
  type = "Combo",
  sample_groups = SCZ_sample_groups
)
# Run SVM with various FD threshold cutoffs using right postcentral cortex catch22
SCZ_ROI_catch22_svm_res <- run_repeat_cv_linear_svm(movement_data = SCZ_movement_data,
                                                    catch22_data = SCZ_catch22_zscored,
                                                    movement_var = "FD",
                                                    type = "Brain Region",
                                                    input_region = "ctx-rh-postcentral",
                                                    sample_groups = SCZ_sample_groups
                                                    
)
# Run SVM with various FD threshold cutoffs using just movement data
SCZ_mvmt_svm_res <- run_repeat_cv_linear_svm(movement_data = SCZ_movement_data,
                                             catch22_data = SCZ_catch22_zscored,
                                             movement_var = "FD",
                                             type = "Movement Only",
                                             sample_groups = SCZ_sample_groups
)

save(SCZ_combo_catch22_svm_res,
     SCZ_ROI_catch22_svm_res,
     SCZ_mvmt_svm_res,
     file=paste0(SCZ_rdata_path, "SCZ_catch22_SVM_with_Movement_Thresholds.Rdata"))
# load(paste0(SCZ_rdata_path, "SCZ_catch22_SVM_with_Movement_Thresholds.Rdata"))

SCZ_combo_catch22_svm_res %>%
  plyr::rbind.fill(., SCZ_ROI_catch22_svm_res) %>%
  plyr::rbind.fill(., SCZ_mvmt_svm_res) %>%
  mutate(Method = case_when(Method=="ctx-rh-postcentral" ~ "Right Postcentral Gyrus",
                            Method=="catch22 Combo" ~ "catch22 Combo",
                            T ~ Method)) %>%
  mutate(Method = factor(Method, levels = c("Movement Only",
                                            "Right Postcentral Gyrus",
                                            "catch22 Combo"))) %>%
  ggplot(data=., mapping=aes(x=movement_threshold)) +
  scale_x_reverse() +
  geom_ribbon(aes(ymin = meanbacc - sdbacc,
                  ymax = meanbacc + sdbacc,
                  fill = Method),
              alpha=0.2) +
  geom_line(aes(y=meanbacc, color = Method)) +
  xlab("FD Maximum Threshold") +
  ylab("Balanced Accuracy (%)")  +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2),
         fill = guide_legend(nrow = 2))
ggsave(paste0(plot_path, "SCZ_FD_threshold_SVM_balanced_accuracy.png"),
       width = 5, height = 3.5, units="in", dpi=300, bg="white")


# Load ASD catch22 z-scored data
noise_proc <- "FC1000"

ASD_catch22_zscored <- readRDS(paste0(ASD_rdata_path, 
                                      "ABIDE_ASD_catch22_filtered_zscored.Rds")) %>%
  filter(Noise_Proc == noise_proc)

# Get diagnosis proportions
ASD_sample_groups <- readRDS(paste0(ASD_rdata_path, 
                                    sprintf("ABIDE_ASD_samples_with_univariate_%s_and_pairwise_%s_filtered.Rds",
                                            univariate_feature_set,
                                            pairwise_feature_set))) %>%
  left_join(., ASD_subject_metadata) %>%
  distinct(Sample_ID, Diagnosis)


# Run SVM with various FD threshold cutoffs using univariate combo catch22
ASD_combo_catch22_svm_res <- run_repeat_cv_linear_svm(
  mvmt_list = seq(0.05, 0.6, by=0.025),
  movement_data = ASD_movement_data,
  catch22_data = ASD_catch22_zscored,
  movement_var = "mean_displacement",
  type = "Combo",
  sample_groups = ASD_sample_groups
)
# Run SVM with various FD threshold cutoffs using right postcentral cortex catch22
ASD_ROI_catch22_svm_res <- run_repeat_cv_linear_svm(
  mvmt_list = seq(0.05, 0.6, by=0.025),
  movement_data = ASD_movement_data,
  catch22_data = ASD_catch22_zscored,
  movement_var = "mean_displacement",
  type = "Brain Region",
  input_region = "Superior Frontal Gyrus",
  sample_groups = ASD_sample_groups
  
)
# Run SVM with various FD threshold cutoffs using just movement data
ASD_mvmt_svm_res <- run_repeat_cv_linear_svm(
  mvmt_list = seq(0.05, 0.6, by=0.025),
  movement_data = ASD_movement_data,
  catch22_data = ASD_catch22_zscored,
  movement_var = "mean_displacement",
  type = "Movement Only",
  sample_groups = ASD_sample_groups
)

save(ASD_combo_catch22_svm_res,
     ASD_ROI_catch22_svm_res,
     ASD_mvmt_svm_res,
     file=paste0(ASD_rdata_path, "ASD_catch22_SVM_with_Movement_Thresholds.Rdata"))
# load(paste0(ASD_rdata_path, "ASD_catch22_SVM_with_Movement_Thresholds.Rdata"))

ASD_combo_catch22_svm_res %>%
  plyr::rbind.fill(., ASD_ROI_catch22_svm_res) %>%
  plyr::rbind.fill(., ASD_mvmt_svm_res) %>%
  mutate(Method = factor(Method, levels = c("Movement Only",
                                            "Superior Frontal Gyrus",
                                            "catch22 Combo"))) %>%
  ggplot(data=., mapping=aes(x=movement_threshold)) +
  scale_x_reverse() +
  geom_ribbon(aes(ymin = meanbacc - sdbacc,
                  ymax = meanbacc + sdbacc,
                  fill = Method),
              alpha=0.2) +
  geom_line(aes(y=meanbacc, color = Method)) +
  xlab("Movement Maximum Threshold") +
  ylab("Balanced Accuracy (%)") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2),
         fill = guide_legend(nrow = 2))
ggsave(paste0(plot_path, "ASD_FD_threshold_SVM_balanced_accuracy.png"),
       width = 5, height = 3.5, units="in", dpi=300, bg="white")