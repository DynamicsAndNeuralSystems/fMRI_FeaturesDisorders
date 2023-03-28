################################################################################
# Load libraries
################################################################################

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/opt/anaconda3/envs/pyspi/bin/python3"
reticulate::use_python(python_to_use)

library(tidyverse)
library(reticulate)
library(icesTAF)
library(cowplot)
library(ggpubr)
library(ggsignif)
library(patchwork)
library(feather)
library(glue)
library(R.matlab)
library(broom)
library(see)
library(colorspace)
library(scales)
theme_set(theme_cowplot())

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

# DIY rlist::list.append
list.append <- function (.data, ...) 
{
  if (is.list(.data)) {
    c(.data, list(...))
  }
  else {
    c(.data, ..., recursive = FALSE)
  }
}


################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
source(paste0(github_dir, "data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/movement_analysis/")
TAF::mkdir(plot_path)
univariate_feature_set <- "catch24"

study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Comparison_Group = c("Schizophrenia", "Bipolar", "ADHD", "ASD"),
                             Group_Title = c("SCZ", "BPD", "ADHD", "ASD"))

data_path <- "~/data/TS_Feature_Manuscript/"
UCLA_CNP_data_path <- "~/data/UCLA_CNP/"
ABIDE_ASD_data_path <- "~/data/ABIDE_ASD/"

# Load data on subjects we actually used
UCLA_CNP_subjects_used <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/processed_data/UCLA_CNP_filtered_sample_info_AROMA_2P_GMR_catch22_pyspi14.feather")) %>%
  pull(Sample_ID)
ABIDE_ASD_subjects_used <- pyarrow_feather$read_feather(glue("{ABIDE_ASD_data_path}/processed_data/ABIDE_ASD_filtered_sample_info_FC1000_catch22_pyspi14.feather")) %>%
  pull(Sample_ID)

# Load subject metadata
UCLA_CNP_sample_metadata <- feather::read_feather(glue("{UCLA_CNP_data_path}/study_metadata/UCLA_CNP_sample_metadata.feather")) %>%
  filter(Sample_ID %in% UCLA_CNP_subjects_used)
ABIDE_ASD_sample_metadata <- feather::read_feather(glue("{ABIDE_ASD_data_path}/study_metadata/ABIDE_ASD_sample_metadata.feather")) %>%
  filter(Sample_ID %in% ABIDE_ASD_subjects_used)

# Load raw feature data
UCLA_CNP_catch24 <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_catch24_filtered.feather")  %>%
  left_join(., UCLA_CNP_sample_metadata)
ABIDE_ASD_catch24 <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/processed_data/ABIDE_ASD_FC1000_catch24_filtered.feather")  %>%
  left_join(., ABIDE_ASD_sample_metadata)
catch24_info <- read.csv(glue("{github_dir}/data_visualisation/catch24_info.csv"))
UCLA_CNP_pyspi14 <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered.feather")  %>%
  group_by(Sample_ID, SPI) %>%
  summarise(mean_across_brain = mean(value, na.rm=T)) %>%
  left_join(., UCLA_CNP_sample_metadata) %>%
  filter(!is.na(Diagnosis))
ABIDE_ASD_pyspi14 <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/processed_data/ABIDE_ASD_FC1000_pyspi14_filtered.feather")  %>%
  group_by(Sample_ID, SPI) %>%
  summarise(mean_across_brain = mean(value, na.rm=T)) %>%
  left_join(., ABIDE_ASD_sample_metadata) %>%
  filter(!is.na(Diagnosis))
pyspi14_info <- read.csv(glue("{github_dir}/data_visualisation/SPI_info.csv"))

# Load mean framewise displacement data
UCLA_CNP_mean_FD <- read.table(glue("{UCLA_CNP_data_path}/movement_data/fmriprep/UCLA_CNP_mFD.txt"), 
                                       sep=",", colClasses = "character")
ABIDE_ASD_mean_FD <- read.table(glue("{ABIDE_ASD_data_path}/movement_data/fmriprep/ABIDE_ASD_mFD.txt"), 
                                     sep=",", colClasses = "character")
# Load full framewise displacement data
UCLA_CNP_full_FD <- as.data.frame(readMat(glue("{UCLA_CNP_data_path}/movement_data/fmriprep/UCLA_CNP_all_FD.mat"))[[1]])
ABIDE_ASD_full_FD <- as.data.frame(readMat(glue("{ABIDE_ASD_data_path}/movement_data/fmriprep/ABIDE_ASD_all_FD.mat"))[[1]])
colnames(UCLA_CNP_full_FD) <- colnames(ABIDE_ASD_full_FD) <- colnames(UCLA_CNP_mean_FD) <- colnames(ABIDE_ASD_mean_FD) <- c("Sample_ID", "Jenkinson", "Power", "VanDijk")

# Un-list full framewise displacement data
UCLA_CNP_full_FD <- UCLA_CNP_full_FD %>%
  unnest(cols=c(Sample_ID, Jenkinson, Power,VanDijk)) %>%
  unnest(cols=c(Sample_ID, Jenkinson, Power, VanDijk)) %>%
  filter(Sample_ID %in% UCLA_CNP_subjects_used) %>%
  group_by(Sample_ID) %>%
  mutate(Frame_Number = row_number())
ABIDE_ASD_full_FD <- ABIDE_ASD_full_FD %>%
  unnest(cols=c(Sample_ID, Jenkinson, Power,VanDijk)) %>%
  unnest(cols=c(Sample_ID, Jenkinson, Power, VanDijk)) %>%
  filter(Sample_ID %in% ABIDE_ASD_subjects_used) %>%
  group_by(Sample_ID) %>%
  mutate(Frame_Number = row_number())
colnames(UCLA_CNP_full_FD) <- colnames(ABIDE_ASD_full_FD) <- c("Sample_ID", "Jenkinson", "Power", "VanDijk", "Frame_Number")
  
# Set mFD columns as numeric
UCLA_CNP_mean_FD <- UCLA_CNP_mean_FD %>%
  mutate_at(c("Jenkinson", "Power", "VanDijk"), function(x) as.numeric(x)) %>%
  left_join(., UCLA_CNP_sample_metadata) %>%
  filter(Sample_ID %in% UCLA_CNP_subjects_used)
ABIDE_ASD_mean_FD <- ABIDE_ASD_mean_FD %>%
  mutate_at(c("Jenkinson", "Power", "VanDijk"), function(x) as.numeric(x)) %>%
  left_join(., ABIDE_ASD_sample_metadata) %>%
  filter(Sample_ID %in% ABIDE_ASD_subjects_used)

################################################################################
# Compare FD-Power distributions between each case-control comparison
################################################################################

control_color <- "#5BB67B"
group_colors <- c("#573DC7", "#D5492A", "#0F9EA9","#C47B2F")

plot_group_vs_control <- function(FD_dataset,
                                  study,
                                  dx,
                                  dx_title,
                                  group_color) {
  p <- FD_dataset %>%
    filter(Study==study,
           Diagnosis %in% c("Control", dx)) %>%
    mutate(Diagnosis = factor(Diagnosis, levels = c(dx, "Control"))) %>%
    ggplot(data=., mapping=aes(x=Diagnosis, y=Power)) +
    geom_violin(aes(fill=Diagnosis)) +
    geom_boxplot(color="black", fill=NA, width=0.1) +
    geom_signif(test = "wilcox.test",
                comparisons = list(c(dx, "Control")), 
                map_signif_level=TRUE) +
    scale_fill_manual(values = c(group_color, control_color)) +
    scale_y_continuous(expand = c(0,0,0.1,0)) +
    scale_x_discrete(labels = c(dx_title, "Control")) +
    ggtitle(dx_title) +
    ylab("mFD-Power") +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          plot.title = element_text(hjust=0.5, size=14))
  return(p)
}

plots <- 1:4 %>%
  purrr::map(~ plot_group_vs_control(plyr::rbind.fill(UCLA_CNP_mean_FD, ABIDE_ASD_mean_FD),
                                     study = study_group_df$Study[.x],
                                     dx = study_group_df$Comparison_Group[.x],
                                     dx_title = study_group_df$Group_Title[.x],
                                     group_color = group_colors[.x]))

wrap_plots(plots, nrow=1)
ggsave(paste0(plot_path, "mFD_Power_by_Group.png"),
       width = 10, height=2.5, units="in", dpi=300)


################################################################################
# Compare correlation of catch24 features with mFD-Power by study
################################################################################

# Find features that correlate most strongly with head motion per study, 
# plot in a barplot
head_motion_univariate_feature_corrs <- UCLA_CNP_catch24 %>%
  plyr::rbind.fill(., ABIDE_ASD_catch24) %>%
  group_by(Sample_ID, names) %>%
  summarise(mean_across_brain = mean(values, na.rm=T)) %>%
  left_join(., plyr::rbind.fill(UCLA_CNP_mean_FD, ABIDE_ASD_mean_FD)) %>%
  left_join(., plyr::rbind.fill(UCLA_CNP_sample_metadata, ABIDE_ASD_sample_metadata)) %>%
  group_by(Study, names) %>%
  nest() %>%
  mutate(
    test = map(data, ~ cor.test(.x$mean_across_brain, .x$Power, method="spearman")), # S3 list-col
    tidied = map(test, tidy)
  ) %>%
  unnest(tidied) %>%
  select(-data, -test) %>%
  group_by(Study) %>%
  mutate(p_Bonferroni = p.adjust(p.value, method="bonferroni")) %>%
  filter(p_Bonferroni < 0.05) %>%
  arrange(desc(abs(estimate)))

head_motion_pairwise_feature_corrs <- UCLA_CNP_pyspi14 %>%
  plyr::rbind.fill(., ABIDE_ASD_pyspi14) %>%
  left_join(., plyr::rbind.fill(UCLA_CNP_mean_FD, ABIDE_ASD_mean_FD)) %>%
  arrange(SPI) %>%
  group_by(Study, SPI) %>%
  nest() %>%
  mutate(
    test = map(data, ~ cor.test(.x$mean_across_brain, .x$Power, method="spearman")), # S3 list-col
    tidied = map(test, tidy)
  ) %>%
  unnest(tidied) %>%
  select(-data, -test) %>%
  group_by(Study) %>%
  mutate(p_Bonferroni = p.adjust(p.value, method="bonferroni")) %>%
  filter(p_Bonferroni < 0.05) %>%
  arrange(desc(abs(estimate))) %>%
  dplyr::rename("names"="SPI")

merged_time_series_info <- catch24_info %>%
  dplyr::rename("names" = "TS_Feature") %>%
  mutate(Feature_Type = "Univariate") %>%
  plyr::rbind.fill(., pyspi14_info %>% 
                     dplyr::rename("names" = "SPI",
                                   "Figure_Name" = "Nickname") %>%
                     mutate(Feature_Type = "Pairwise"))
# Annotation bar with feature type
plyr::rbind.fill(head_motion_univariate_feature_corrs,
                 head_motion_pairwise_feature_corrs) %>%
  left_join(., merged_time_series_info) %>%
  mutate(Study = ifelse(Study == "ABIDE_ASD", "ABIDE", "UCLA CNP")) %>%
  mutate(Figure_Name = fct_reorder(Figure_Name, estimate, .fun="mean", .desc=F),
         Feature_Type = fct_reorder(Feature_Type, estimate, .fun="mean", .desc=F)) %>%
  ggplot(data=., mapping=aes(x=0, y=Figure_Name, fill=Feature_Type)) +
  geom_tile() +
  theme_void() +
  scale_fill_manual(values=c("#803556", "#E8A6A9")) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text=element_text(size=14)) +
  guides(fill = guide_legend(title.position = "top", 
                             ncol = 2,
                             byrow=T,
                             title.hjust = 0.5)) 
ggsave(glue("{plot_path}/Feature_type_colorbar.png"),
       width=6, height=6, units="in", dpi=300)

# Heatmap
plyr::rbind.fill(head_motion_univariate_feature_corrs,
                 head_motion_pairwise_feature_corrs) %>%
  left_join(., merged_time_series_info) %>%
  mutate(Study = ifelse(Study == "ABIDE_ASD", "ABIDE", "UCLA CNP")) %>%
  mutate(Figure_Name = fct_reorder(Figure_Name, estimate, .fun="mean")) %>%
  ggplot(data=., mapping=aes(x=Study, y=Figure_Name, fill=estimate)) +
  geom_tile() +
  geom_text(aes(label = round(estimate,2)),
            size=5) +
  scale_fill_continuous_divergingx(palette="RdBu", rev=TRUE) +
  ylab("Time-series feature") +
  labs(fill="mFD \u03C1")
ggsave(glue("{plot_path}/Movement_spearman_feature_corr.png"),
       width=6, height=6, units="in", dpi=300)


################################################################################
# Plot individual features vs movement
plot_feature_vs_movement <- function(feature_name, dataset_to_use, 
                                     title_label, y_label, rho_pos) {
  feature_data <- UCLA_CNP_catch24 %>%
    plyr::rbind.fill(., ABIDE_ASD_catch24) %>%
    filter(names == feature_name) %>%
    group_by(Sample_ID) %>%
    summarise(meanval = mean(values, na.rm=T)) %>%
    left_join(., plyr::rbind.fill(UCLA_CNP_mean_FD, ABIDE_ASD_mean_FD)) %>%
    left_join(., plyr::rbind.fill(UCLA_CNP_sample_metadata, ABIDE_ASD_sample_metadata)) %>%
    mutate(Study = ifelse(Study == "UCLA_CNP", "UCLA CNP", "ABIDE"))  %>%
    mutate(Diagnosis = case_when(Diagnosis == "Schizophrenia" ~ "SCZ",
                                        Diagnosis == "Bipolar" ~ "BPD",
                                        T ~ Diagnosis))
  
  p <- feature_data %>%
    ggplot(data=., mapping=aes(x=Power, y=meanval)) +
    xlab("mFD-Power") +
    ylab(y_label) +
    ggtitle(title_label) +
    facet_grid(Study ~ ., scales="free", switch="both") +
    geom_point(aes(color=Diagnosis)) +
    scale_color_manual(values=c("Control" = "#5BB67B", 
                                "SCZ" = "#573DC7", 
                                "BPD" = "#D5492A", 
                                "ADHD" = "#0F9EA9", 
                                "ASD" = "#C47B2F")) +
    stat_smooth(method="lm", 
                data=subset(feature_data, Study == dataset_to_use),
                color="black",
                se=F) +
    stat_cor(method="spearman", 
             aes(label = after_stat(r.label)),
             data=subset(feature_data, Study == dataset_to_use),
             size = 6,
             label.x = rho_pos,
             cor.coef.name = "rho", 
             p.accuracy = 0) +
    guides(color = guide_legend(title.position = "top", 
                                nrow = 3,
                                byrow=T,
                                override.aes = list(size=5),
                                title.hjust = 0.5)) +
    theme(legend.position="bottom",
          plot.title = element_text(hjust=0.5),
          strip.placement = "outside")
  
  return(p)
}

# SD
plot_feature_vs_movement(feature_name = "DN_Spread_Std",
                         dataset_to_use = "ABIDE",
                         title_label = "Standard deviation",
                         y_label = "Mean SD across brain",
                         rho_pos = 0)
ggsave(paste0(plot_path, "mFD_vs_SD.png"),
       width = 3.75, height=5.75, units="in", dpi=300)

# CO_Embed2_Dist_tau_d_expfit_meandiff
plot_feature_vs_movement(feature_name = "CO_Embed2_Dist_tau_d_expfit_meandiff",
                         dataset_to_use = "UCLA CNP",
                         title_label = "Embedding distance",
                         y_label = "Mean embedding distance across brain",
                         rho_pos = 0.35)
ggsave(paste0(plot_path, "mFD_vs_embedding_distance.png"),
       width = 3.75, height=5.75, units="in", dpi=300)


################################################################################
# Compare balanced accuracy using just movement data for SVM
################################################################################

# Load null data for comparison
univariate_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_null_balanced_accuracy_distributions.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set,
         Analysis_Type == "Univariate_Combo") %>%
  dplyr::select(Study, Comparison_Group, Null_Balanced_Accuracy) %>%
  dplyr::rename("Repeat_Balanced_Accuracy" = "Null_Balanced_Accuracy") %>%
  mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BPD",
                                      T ~ Comparison_Group)) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD")))

# Import normalisation methods from python
source_python("~/github/fMRI_FeaturesDisorders/helper_functions/classification/mixed_sigmoid_normalisation.py")

# Import python SVM functions
np <- import("numpy", convert=FALSE)
sklearn_SVM <- import("sklearn.svm")
sklearn_pipeline <- import("sklearn.pipeline")
sklearn_modsel <- import("sklearn.model_selection")

# Define standard model parameters
pipe <- sklearn_pipeline$Pipeline(list(c('scaler', MixedSigmoidScaler(unit_variance=TRUE)), 
                                       c('SVM', 
                                         sklearn_SVM$SVC(kernel = "linear", C = 1, 
                                                         shrinking = FALSE, 
                                                         class_weight = "balanced"))))

# Function to run 10-repeat 10-fold CV for each comparison group
run_movement_SVM <- function(movement_dataset,
                             comparison_group,
                             pipe_obj) {
  # Comparison group vs. Controls
  group_ctrl_mvmt_data <- movement_dataset %>%
    filter(Diagnosis %in% c("Control", comparison_group)) %>%
    dplyr::select(Sample_ID, Diagnosis, Power) %>%
    mutate(Power = as.numeric(Power)) %>%
    as.matrix()
  
  # Run 10 repeats of 10-fold CV using just movement data
  repeat_movement_balacc_list <- list()
  
  for (repeat_i in 1:10) {
    skf = sklearn_modsel$StratifiedKFold(n_splits=as.integer(10), shuffle=TRUE, random_state=as.integer(repeat_i))
    movement_cv_results_for_repeat <- sklearn_modsel$cross_validate(pipe_obj, 
                                                                    np$array(as.numeric(group_ctrl_mvmt_data[,3]))$reshape(as.integer(-1),
                                                                                                                           as.integer(1)), 
                                                                    np$array(group_ctrl_mvmt_data[,2]), 
                                                                    cv=skf, 
                                                                    scoring="balanced_accuracy",
                                                                    n_jobs = 1,
                                                                    return_estimator=FALSE)$test_score
    
    # Compile results
    repeat_res <- data.frame(Comparison_Group = comparison_group,
                             Repeat = repeat_i,
                             Balanced_Accuracy = movement_cv_results_for_repeat,
                             Fold = 1:10)
    # Append to list
    repeat_movement_balacc_list <- list.append(repeat_movement_balacc_list, repeat_res)
    
  }
  # Integreate results across repeats
  movement_results_for_group <- do.call(plyr::rbind.fill, repeat_movement_balacc_list)
  
}

# Schizophrenia
scz_mvmt_res <- run_movement_SVM(movement_dataset = UCLA_CNP_mean_FD,
                                 comparison_group = "Schizophrenia",
                                 pipe_obj = pipe)
# Bipolar
bpd_mvmt_res <- run_movement_SVM(movement_dataset = UCLA_CNP_mean_FD,
                                 comparison_group = "Bipolar",
                                 pipe_obj = pipe)

# ADHD
adhd_mvmt_res <- run_movement_SVM(movement_dataset = UCLA_CNP_mean_FD,
                                  comparison_group = "ADHD",
                                  pipe_obj = pipe)

# ASD
asd_mvmt_res <- run_movement_SVM(movement_dataset = ABIDE_ASD_mean_FD,
                                 comparison_group = "ASD",
                                 pipe_obj = pipe)

# Combine results
do.call(plyr::rbind.fill, list(scz_mvmt_res, bpd_mvmt_res, adhd_mvmt_res, asd_mvmt_res)) %>%
  group_by(Comparison_Group, Repeat) %>%
  summarise(Repeat_Balanced_Accuracy = mean(Balanced_Accuracy, na.rm=T)) %>%
  group_by(Comparison_Group) %>%
  mutate(Balanced_Accuracy_Across_Repeats = mean(Repeat_Balanced_Accuracy, na.rm=T),
         Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BPD",
                                      T ~ Comparison_Group)) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD")))%>%
  ggplot(data=., mapping=aes(x=Comparison_Group, y=100*Repeat_Balanced_Accuracy)) +
  geom_violinhalf(aes(fill=Comparison_Group), scale="width")  +
  geom_boxplot(width=0.1, notch=FALSE, notchwidth = 0.4, outlier.shape = NA,
               fill=NA, color=alpha("white", 0.7),
               position = position_nudge(x=0.058), coef = 0) +
  geom_violinhalf(data = univariate_null_distribution, fill="gray80", flip=T) +
  geom_boxplot(data = univariate_null_distribution, 
               width=0.1, notch=FALSE, 
               notchwidth = 0.4, outlier.shape = NA,
               fill=NA, color=alpha("black", 0.7),
               position = position_nudge(x=-0.058), coef = 0) +
  geom_hline(yintercept = 50, linetype=2, alpha=0.5) +
  scale_fill_manual(values=c("Control" = "#5BB67B", 
                              "SCZ" = "#573DC7", 
                              "BPD" = "#D5492A", 
                              "ADHD" = "#0F9EA9", 
                              "ASD" = "#C47B2F")) +
  ylab("Mean Balanced Accuracy\nby Repeat (%)") +
  xlab("Clinical Group") +
  theme(legend.position = "none")
ggsave(paste0(plot_path, "SVM_movement_balacc.png"),
       width = 3, height=3, units="in", dpi=300)
