################################################################################
# Load libraries
################################################################################

library(tidyverse)
library(icesTAF)
library(cowplot)
library(ggpubr)
library(ggsignif)
library(patchwork)
library(feather)
library(glue)
theme_set(theme_cowplot())

################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
source(paste0(github_dir, "data_visualisation/manuscript_figures/Manuscript_Draft_Visualisations_Helper.R"))
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/FigureS1/")
TAF::mkdir(plot_path)

UCLA_CNP_data_path <- "~/data/UCLA_CNP/"
ABIDE_ASD_data_path <- "~/data/ABIDE_ASD/"

# Load subject metadata
UCLA_CNP_sample_metadata <- feather::read_feather(glue("{UCLA_CNP_data_path}/study_metadata/UCLA_CNP_sample_metadata.feather"))
ABIDE_ASD_sample_metadata <- feather::read_feather(glue("{ABIDE_ASD_data_path}/study_metadata/ABIDE_ASD_sample_metadata.feather"))

# Load movement data
UCLA_CNP_movement_data <- read.table(glue("{UCLA_CNP_data_path}/movement_data/UCLA_CNP_mFD.txt"), 
                                       sep=",", colClasses = "character")
ABIDE_ASD_movement_data <- read.table(glue("{ABIDE_ASD_data_path}/movement_data/ABIDE_ASD_mFD.txt"), 
                                     sep=",", colClasses = "character")
colnames(UCLA_CNP_movement_data) <- colnames(ABIDE_ASD_movement_data) <- c("Sample_ID", "Jenkinson", "Power", "VanDijk")

# Set mFD columns as numeric
UCLA_CNP_movement_data <- UCLA_CNP_movement_data %>%
  mutate_at(c("Jenkinson", "Power", "VanDijk"), function(x) as.numeric(x)) %>%
  left_join(., UCLA_CNP_sample_metadata)
ABIDE_ASD_movement_data <- ABIDE_ASD_movement_data %>%
  mutate_at(c("Jenkinson", "Power", "VanDijk"), function(x) as.numeric(x)) %>%
  left_join(., ABIDE_ASD_sample_metadata)

################################################################################
# Compare FD-Power distributions between each case-control comparison
################################################################################'

# Is age related to movement in either cohort?
UCLA_CNP_movement_data %>%
  mutate(Cohort = "UCLA CNP", Age = as.numeric(Age)) %>%
  ggplot(data=., mapping=aes(x=Age, y=Power, color=Diagnosis)) +
  geom_point() +
  geom_smooth(method="lm")

ABIDE_ASD_movement_data %>%
  mutate(Cohort = "ABIDE ASD", Age = as.numeric(Age)) %>%
  filter(Power< 3) %>%
  ggplot(data=., mapping=aes(x=Age, y=Power, color=Diagnosis)) +
  geom_point() +
  geom_smooth(method="lm")


# Maybe we only focus on the FDpower estimates for the paper?
UCLA_CNP_movement_data %>%
  mutate(Cohort = "UCLA CNP") %>%
  ggplot(data=., mapping=aes(x=Diagnosis, y=Power)) +
  geom_violin(aes(fill=Diagnosis)) +
  geom_boxplot(color="black", fill=NA, width=0.1) +
  facet_wrap(Cohort ~ ., scales="free_x") +
  geom_signif(test = "wilcox.test",
              comparisons = list(c("Schizophrenia", "Control")), 
              map_signif_level=TRUE) +
  scale_fill_manual(values = c("#00B06D", "#737373")) +
  scale_y_continuous(expand = c(0,0,0.15,0)) +
  ylab("Head Movement\n(mFD-Power)") +
  xlab("Group") +
  theme(legend.position = "none")
ggsave(paste0(plot_path, "UCLA_CNP_mFD_Power_by_Group.png"),
       width = 3, height=2.25, units="in", dpi=300, bg="white")

ABIDE_ASD_movement_data %>%
  filter(Power > 1)
ABIDE_ASD_movement_data %>%
  filter(Power < 1) %>%
  mutate(Cohort = "ABIDE") %>%
  ggplot(data=., mapping=aes(x=Diagnosis, y=Power)) +
  geom_violin(aes(fill=Diagnosis)) +
  geom_boxplot(color="black", fill=NA, width=0.1) +
  facet_wrap(Cohort ~ ., scales="free_x") +
  geom_signif(test = "wilcox.test",
              comparisons = list(c("ASD", "Control")), 
              map_signif_level=TRUE) +
  scale_fill_manual(values = c("#00B06D", "#737373")) +
  scale_y_continuous(expand = c(0,0,0.15,0)) +
  ylab("Head Movement\n(mFD-Power)") +
  xlab("Group") +
  theme(legend.position = "none")
ggsave(paste0(plot_path, "ABIDE_mFD_Power_by_Group.png"),
       width = 3, height=2.25, units="in", dpi=300, bg="white")


################################################################################
# Plot # subjects retained per group by motion threshold
################################################################################
# Function to find number of subjects retained per group by motion threshold
subjects_retained_by_motion <- function(movement_data, 
                                        motion_variable = "Power",
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
UCLA_subjects_by_motion <- subjects_retained_by_motion(movement_data = subset(UCLA_ABIDE_movement_data, Study=="UCLA_CNP"),
                                                      motion_variable = "Power",
                                                      motion_range=seq(0, max(subset(UCLA_ABIDE_movement_data, Study=="UCLA_CNP") %>%
                                                                                pull(Power)), by=0.02))
ABIDE_subjects_by_motion <- subjects_retained_by_motion(movement_data = subset(UCLA_ABIDE_movement_data, Study=="ABIDE_ASD"),
                                                      motion_variable = "Power",
                                                      motion_range=seq(0, max(subset(UCLA_ABIDE_movement_data,
                                                                                     Study=="ABIDE_ASD") %>%
                                                                                pull(Power)), by=0.02))

# Find threshold above which 75+ % of subjects are retained per group
UCLA_range_to_evaluate <- UCLA_subjects_by_motion %>%
  group_by(Diagnosis) %>%
  # Find number of subjects needed to maintain 75% of sample size
  mutate(Total_N_by_Group = max(Num_Subjects),
         Threshold = ceiling(0.75*Total_N_by_Group)) %>% 
  ungroup() %>%
  group_by(Motion_Threshold) %>%
  filter(all(Num_Subjects >= Threshold)) %>%
  distinct(Motion_Threshold) %>%
  pull(Motion_Threshold)

ABIDE_range_to_evaluate <- ABIDE_subjects_by_motion %>%
  group_by(Diagnosis) %>%
  # Find number of subjects needed to maintain 75% of sample size
  mutate(Total_N_by_Group = max(Num_Subjects),
         Threshold = ceiling(0.75*Total_N_by_Group)) %>% 
  ungroup() %>%
  group_by(Motion_Threshold) %>%
  filter(all(Num_Subjects >= Threshold)) %>%
  distinct(Motion_Threshold) %>%
  # Only look up to an FD power of 1
  filter(Motion_Threshold <= 1) %>%
  pull(Motion_Threshold)

# Plot # of subjects retained per FD threshold by group
# SCZ
UCLA_num_motion <- UCLA_subjects_by_motion %>%
  group_by(Diagnosis) %>%
  ggplot(data=., mapping=aes(x=Motion_Threshold, y=Num_Subjects, 
                             color=Diagnosis, group=Diagnosis)) +
  geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=min(UCLA_range_to_evaluate), 
                xmax=max(UCLA_range_to_evaluate)),
            color=NA, fill="gray85", alpha=0.5) +
  geom_line(size=1.75, alpha=0.9) +
  ylab("# Subjects") +
  scale_color_manual(values = c("#00B06D", "#737373")) +
  scale_x_reverse(limits=c(1,0)) +
  theme(legend.position="none",
        axis.title.x = element_blank(),
        strip.placement = "outside",
        plot.title=element_text(hjust=0.5))
# ASD
ABIDE_num_motion <- ABIDE_subjects_by_motion %>%
  mutate(Diagnosis = factor(Diagnosis, levels=c("Control", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Motion_Threshold, y=Num_Subjects, 
                             color=Diagnosis, group=Diagnosis)) +
  geom_rect(aes(ymin=-Inf, ymax=Inf, 
                xmin=min(ABIDE_range_to_evaluate),
                xmax=max(ABIDE_range_to_evaluate)),
            color=NA, fill="gray85", alpha=0.5) +
  geom_line(linewidth=1.75, alpha=0.8) +
  ylab("# Subjects") +
  scale_color_manual(values = c("#00B06D", "#737373")) +
  scale_x_reverse(limits=c(1,0)) +
  theme(legend.position="none",
        axis.title.x = element_blank(),
        strip.placement = "outside",
        plot.title=element_text(hjust=0.5))

# Plot % of subjects retained per FD threshold by group
# SCZ
UCLA_perc_motion <- UCLA_subjects_by_motion %>%
  group_by(Diagnosis) %>%
  mutate(Perc_Subjects = Num_Subjects / max(Num_Subjects, na.rm=T)) %>%
  ggplot(data=., mapping=aes(x=Motion_Threshold, y=100*Perc_Subjects, 
                             color=Diagnosis, group=Diagnosis)) +
  geom_rect(aes(ymin=-Inf, ymax=Inf, 
                xmin=min(UCLA_range_to_evaluate), 
                xmax=max(UCLA_range_to_evaluate)),
            color=NA, fill="gray85", alpha=0.5) +
  geom_line(size=1.75, alpha=0.9) +
  ylab("% of Subjects") +
  labs(color="Group") +
  xlab("Movement (mFD-Power)\nMaximum Threshold") +
  scale_color_manual(values = c("#00B06D", "#737373")) +
  scale_x_reverse(limits=c(1,0)) +
  theme(legend.position="bottom",
        strip.placement = "outside",
        plot.title=element_text(hjust=0.5))
# ASD
ABIDE_perc_motion <- ABIDE_subjects_by_motion %>%
  mutate(Diagnosis = factor(Diagnosis, levels=c("Control", "ASD"))) %>%
  group_by(Diagnosis) %>%
  mutate(Perc_Subjects = Num_Subjects / max(Num_Subjects, na.rm=T)) %>%
  ggplot(data=., mapping=aes(x=Motion_Threshold, y=100*Perc_Subjects, 
                             color=Diagnosis, group=Diagnosis)) +
  geom_rect(aes(ymin=-Inf, ymax=Inf, 
                xmin=min(ABIDE_range_to_evaluate),
                xmax=max(ABIDE_range_to_evaluate)),
            color=NA, fill="gray85", alpha=0.5) +
  geom_line(size=1.75, alpha=0.9) +
  ylab("% of Subjects") +
  labs(color="Group") +
  xlab("Movement (mFD-Power)\nMaximum Threshold") +
  scale_color_manual(values = c("#00B06D", "#737373")) +
  scale_x_reverse(limits=c(1,0)) +
  theme(legend.position="bottom",
        strip.placement = "outside",
        plot.title=element_text(hjust=0.5))

ABIDE_num_motion/ABIDE_perc_motion
ggsave(paste0(plot_path, "ABIDE_mvmt_threshold_num_subjects.png"),
       width = 4, height = 4, units="in", dpi=300)

UCLA_num_motion/UCLA_perc_motion
ggsave(paste0(plot_path, "UCLA_mvmt_threshold_num_subjects.png"),
       width = 4, height = 4, units="in", dpi=300)

################################################################################
# Movement-based SVM results
################################################################################

UCLA_combo_catch22_svm_res %>%
  plyr::rbind.fill(., UCLA_ROI_catch22_svm_res) %>%
  plyr::rbind.fill(., UCLA_mvmt_svm_res) %>%
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
ggsave(paste0(plot_path, "UCLA_FD_threshold_SVM_balanced_accuracy.png"),
       width = 5, height = 3.5, units="in", dpi=300, bg="white")


ABIDE_combo_catch22_svm_res %>%
  plyr::rbind.fill(., ABIDE_ROI_catch22_svm_res) %>%
  plyr::rbind.fill(., ABIDE_mvmt_svm_res) %>%
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
ggsave(paste0(plot_path, "ABIDE_FD_threshold_SVM_balanced_accuracy.png"),
       width = 5, height = 3.5, units="in", dpi=300, bg="white")