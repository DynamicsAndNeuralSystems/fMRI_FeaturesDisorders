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

data_path <- "~/data/UCLA_CNP_ABIDE_ASD/"
rdata_path <- paste0(data_path, "processed_data/Rdata/")

# Load subject metadata
sample_metadata <- readRDS(paste0(data_path, "study_metadata/UCLA_CNP_ABIDE_ASD_sample_metadata.Rds"))

# Load movement data
UCLA_ABIDE_movement_data <- readRDS(paste0(data_path, "movement_data/UCLA_CNP_ABIDE_ASD_FD.Rds"))
UCLA_movement_data_Linden <- readRDS(paste0(data_path, "movement_data/UCLA_CNP_Linden_FD_Power.Rds"))

################################################################################
# Compare distributions of each mFD type by dataset
################################################################################

# NOTE: I'm plotting the 1 to the 95th percentiles for each dataset
merged_midrange_data <- UCLA_ABIDE_movement_data %>%
  pivot_longer(cols=c(Jenkinson:VanDijk),
               names_to = "Method",
               values_to = "mFD") %>%
  group_by(Cohort, Method) %>%
  mutate(perc025 = quantile(mFD, probs=c(0.025))) %>%
  mutate(perc975 = quantile(mFD, probs=c(0.975))) %>%
  ungroup() %>%
  rowwise() %>%
  filter(perc025 <= mFD, mFD <= perc975) %>%
  mutate(Diagnosis = factor(Diagnosis, levels = c("Control", "ASD", "Schizophrenia")))

merged_midrange_data %>%
  ggplot(data=., mapping=aes(x=Diagnosis, y=mFD)) +
  geom_violin(aes(fill=Diagnosis)) +
  geom_boxplot(color="black", fill=NA, width=0.1) +
  facet_wrap(Cohort ~ Method, scales="free") +
  geom_signif(data = subset(merged_midrange_data, Cohort=="ASD Study"),
              test = "wilcox.test",
              comparisons = list(c("ASD", "Control")), 
              map_signif_level=TRUE) +
  geom_signif(data = subset(merged_midrange_data, Cohort=="SCZ Study"),
              test = "wilcox.test",
              comparisons = list(c("Schizophrenia", "Control")), 
              map_signif_level=TRUE) +
  scale_y_continuous(expand = c(0,0,0.15,0)) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank())
ggsave(paste0(plot_path, "UCLA_ABIDE_mFD_by_Group.png"),
       width = 6, height=4, units="in", dpi=300, bg="white")

# Compare control subjects per cohort
merged_midrange_data %>%
  filter(Diagnosis=="Control") %>%
  ggplot(data=., mapping=aes(x=mFD)) +
  stat_density(aes(fill = Cohort, y=after_stat(scaled)),
               position = "identity", alpha=0.6) +
  ggtitle("mFD Values in Control Participants") +
  ylab("Scaled Density\nof Control Participants") +
  xlab("mFD By Method") +
  facet_wrap(. ~ Method, scales="free") +
  scale_x_continuous(breaks = scales::pretty_breaks(3)) +
  theme(plot.title = element_text(hjust=0.5),
        legend.position="bottom")
ggsave(paste0(plot_path, "Control_mFD_Distribution.png"),
       width = 7, height=3, units="in", dpi=300, bg="white")

# Maybe we only focus on the FDpower estimates for the paper?
ABIDE_power <- merged_midrange_data %>%
  filter(Method == "Power",
         Cohort == "ABIDE Study") 

UCLA_power <- merged_midrange_data %>%
  filter(Method == "Power",
         Cohort == "UCLA Study") 

ABIDE_power %>%
  ggplot(data=., mapping=aes(x=Diagnosis, y=mFD)) +
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

UCLA_power  %>%
  ggplot(data=., mapping=aes(x=Diagnosis, y=mFD)) +
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
ggsave(paste0(plot_path, "UCLA_mFD_Power_by_Group.png"),
       width = 3, height=2.25, units="in", dpi=300, bg="white")

################################################################################
# Compare Linden's mFD-Jenkinson values with Annie's
################################################################################

UCLA_movement_data_Linden %>%
  ggplot(data=., mapping=aes(x=Jenkinson_Linden, y=Power)) +
  geom_point() +
  geom_abline(color="red", slope=1, intercept=0, alpha=0.6, linewidth=1) +
  coord_equal() +
  ggtitle("Mean Framewise Displacement with\nPower et al. 2012 Method") +
  ylab("mFD-Power, Annie") +
  xlab("mFD-Power, Linden") +
  theme(plot.title=element_text(hjust=0.5, size=12))
ggsave(paste0(plot_path, "UCLA_mFD_Power_Linden_Annie.png"),
       width = 4, height=4, units="in", dpi=300, bg="white")

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
                                                      motion_range=seq(0, 1.25, by=0.02))
ABIDE_subjects_by_motion <- subjects_retained_by_motion(movement_data = subset(UCLA_ABIDE_movement_data, Study=="ABIDE_ASD"),
                                                      motion_variable = "Power",
                                                      motion_range=seq(0, 1.25, by=0.02))


# Plot # of subjects retained per FD threshold by group
# SCZ
UCLA_num_motion <- UCLA_subjects_by_motion %>%
  group_by(Diagnosis) %>%
  ggplot(data=., mapping=aes(x=Motion_Threshold, y=Num_Subjects, 
                             color=Diagnosis, group=Diagnosis)) +
  geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=0.12, xmax=0.5),
            color=NA, fill="gray85", alpha=0.5) +
  geom_line(size=1.75, alpha=0.9) +
  ylab("# Subjects") +
  scale_color_manual(values = c("#00B06D", "#737373")) +
  scale_x_reverse(limits=c(1.2,0)) +
  theme(legend.position="none",
        axis.title.x = element_blank(),
        strip.placement = "outside",
        plot.title=element_text(hjust=0.5))
# ASD
ABIDE_num_motion <- ABIDE_subjects_by_motion %>%
  mutate(Diagnosis = factor(Diagnosis, levels=c("Control", "ASD"))) %>%
  ggplot(data=., mapping=aes(x=Motion_Threshold, y=Num_Subjects, 
                             color=Diagnosis, group=Diagnosis)) +
  geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=0.01, xmax=1),
            color=NA, fill="gray85", alpha=0.5) +
  geom_line(size=1.75, alpha=0.8) +
  ylab("# Subjects") +
  scale_color_manual(values = c("#00B06D", "#737373")) +
  scale_x_reverse(limits=c(1.2,0)) +
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
  geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=0.12, xmax=0.6),
            color=NA, fill="gray85", alpha=0.5) +
  geom_line(size=1.75, alpha=0.9) +
  ylab("% of Subjects") +
  labs(color="Group") +
  xlab("Movement (mFD-Power)\nMaximum Threshold") +
  scale_color_manual(values = c("#00B06D", "#737373")) +
  scale_x_reverse(limits=c(1.2,0)) +
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
  geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=0.01, xmax=1),
            color=NA, fill="gray85", alpha=0.5) +
  geom_line(size=1.75, alpha=0.9) +
  ylab("% of Subjects") +
  labs(color="Group") +
  xlab("Movement (mFD-Power)\nMaximum Threshold") +
  scale_color_manual(values = c("#00B06D", "#737373")) +
  scale_x_reverse(limits=c(1.2,0)) +
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