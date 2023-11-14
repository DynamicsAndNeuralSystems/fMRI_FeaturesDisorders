################################################################################
# Load libraries
################################################################################

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"
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
library(ggseg)
library(ggsegHO)
library(ggsegDefaultExtra)
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

github_dir <- "~/Library/CloudStorage/OneDrive-TheUniversityofSydney(Students)/github/fMRI_FeaturesDisorders/"

source(paste0(github_dir, "data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/movement_analysis/")
TAF::mkdir(plot_path)
univariate_feature_set <- "catch25"

study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Noise_Proc = c(rep("AROMA+2P+GMR",3), "FC1000"),
                             Comparison_Group = c("Schizophrenia", "Bipolar", "ADHD", "ASD"),
                             Group_Nickname = c("SCZ", "BP", "ADHD", "ASD"))

data_path <- "~/data/TS_feature_manuscript/"
UCLA_CNP_data_path <- "~/data/UCLA_CNP/"
ABIDE_ASD_data_path <- "~/data/ABIDE_ASD/"

# Load in feature info
catch25_info <- read.csv(glue("{github_dir}/data_visualisation/catch25_info.csv"))
pyspi14_info <- read.csv(glue("{github_dir}/data_visualisation/SPI_info.csv"))

# Load data on subjects we actually used
UCLA_CNP_subjects_used <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/processed_data/UCLA_CNP_filtered_sample_info_AROMA_2P_GMR_catch25_pyspi14.feather")) %>%
  pull(Sample_ID)
ABIDE_ASD_subjects_used <- pyarrow_feather$read_feather(glue("{ABIDE_ASD_data_path}/processed_data/ABIDE_ASD_filtered_sample_info_FC1000_catch25_pyspi14.feather")) %>%
  pull(Sample_ID)

# Load subject metadata
UCLA_CNP_sample_metadata <- feather::read_feather(glue("{UCLA_CNP_data_path}/study_metadata/UCLA_CNP_sample_metadata.feather")) %>%
  filter(Sample_ID %in% UCLA_CNP_subjects_used)
ABIDE_ASD_sample_metadata <- feather::read_feather(glue("{ABIDE_ASD_data_path}/study_metadata/ABIDE_ASD_sample_metadata.feather")) %>%
  filter(Sample_ID %in% ABIDE_ASD_subjects_used)

# Load mean framewise displacement data
UCLA_CNP_mean_FD <- read.table(glue("{UCLA_CNP_data_path}/movement_data/fmriprep/UCLA_CNP_mFD.txt"), 
                                       sep=",", colClasses = "character")
ABIDE_ASD_mean_FD <- read.table(glue("{ABIDE_ASD_data_path}/movement_data/fmriprep/ABIDE_ASD_mFD.txt"),
                                     sep=",", colClasses = "character")
colnames(UCLA_CNP_mean_FD) <- colnames(ABIDE_ASD_mean_FD) <- c("Sample_ID", "Jenkinson", "Power", "VanDijk")

# Set mFD columns as numeric
UCLA_CNP_mean_FD <- UCLA_CNP_mean_FD %>%
  mutate_at(c("Jenkinson", "Power", "VanDijk"), function(x) as.numeric(x)) %>%
  left_join(., UCLA_CNP_sample_metadata, by="Sample_ID") %>%
  filter(Sample_ID %in% UCLA_CNP_subjects_used) %>%
  mutate(Study = "UCLA_CNP")

ABIDE_ASD_mean_FD <- ABIDE_ASD_mean_FD %>%
  mutate_at(c("Jenkinson", "Power", "VanDijk"), function(x) as.numeric(x)) %>%
  left_join(., ABIDE_ASD_sample_metadata) %>%
  filter(Sample_ID %in% ABIDE_ASD_subjects_used) %>%
  mutate(Study = "ABIDE_ASD")

# Load catch25 data
UCLA_CNP_catch25 <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_catch25_filtered.feather")  %>%
  left_join(., UCLA_CNP_sample_metadata)
ABIDE_ASD_catch25 <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/processed_data/ABIDE_ASD_FC1000_catch25_filtered.feather")  %>%
  left_join(., ABIDE_ASD_sample_metadata)

################################################################################
# Compare FD-Power distributions between each case-control comparison
################################################################################

control_color <- "#439E47"
group_colors <- c("#9d60a8", "#2F77C0", "#e45075","#E28328")

plot_group_vs_control <- function(FD_dataset,
                                  study,
                                  dx,
                                  dx_title,
                                  ymin,
                                  ymax,
                                  group_color) {
  p <- FD_dataset %>%
    filter(Study==study,
           Diagnosis %in% c("Control", dx)) %>%
    mutate(Diagnosis = factor(Diagnosis, levels = c(dx, "Control"))) %>%
    ggplot(data=., mapping=aes(x=Diagnosis, y=as.numeric(Power))) +
    geom_violinhalf(aes(fill=Diagnosis), scale="width", 
                    position = position_nudge(x=0.2))  +
    geom_boxplot(width=0.1, notch=FALSE, notchwidth = 0.4, outlier.shape = NA,
                 fill=NA, color="white",
                 position = position_nudge(x=0.27), coef = 0) +
    geom_point(aes(color = Diagnosis), position = position_jitterdodge(dodge.width = 1,
                                                                       jitter.width = 0.5),
               size = 1, alpha=0.7) +
    geom_signif(test = "wilcox.test",
                comparisons = list(c(dx, "Control")), 
                y_position = c(0.55, 0.6),
                map_signif_level=TRUE) +
    scale_fill_manual(values = c(group_color, control_color)) +
    scale_color_manual(values = c(group_color, control_color)) +
    scale_y_continuous(limits = c(ymin, ymax), expand = c(0,0,0.1,0)) +
    scale_x_discrete(labels = c(dx_title, "Control")) +
    ggtitle(dx_title) +
    ylab("mFD-Power") +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          plot.title = element_text(hjust=0.5, size=14))
  return(p)
}

# Find ymin and ymax
FD_dataset = plyr::rbind.fill(UCLA_CNP_mean_FD, ABIDE_ASD_mean_FD)
ymin <- 0
ymax <- max(FD_dataset$Power)
plots <- 1:4 %>%
  purrr::map(~ plot_group_vs_control(FD_dataset = FD_dataset,
                                     study = study_group_df$Study[.x],
                                     dx = study_group_df$Comparison_Group[.x],
                                     dx_title = study_group_df$Group_Nickname[.x],
                                     ymin = ymin,
                                     ymax = 0.61,
                                     group_color = group_colors[.x]))

wrap_plots(plots, nrow=1)
ggsave(paste0(plot_path, "mFD_Power_by_Group.svg"),
       width = 10, height=2.5, units="in", dpi=300)

################################################################################
# Compare region-specific SD values with whole-brain head movement
################################################################################

# Calculate correlations
head_motion_SD_corr_regional <- UCLA_CNP_catch25 %>%
  plyr::rbind.fill(., ABIDE_ASD_catch25) %>%
  filter(names == "DN_Spread_Std") %>%
  left_join(., plyr::rbind.fill(UCLA_CNP_mean_FD, ABIDE_ASD_mean_FD)) %>%
  left_join(., plyr::rbind.fill(UCLA_CNP_sample_metadata, ABIDE_ASD_sample_metadata)) %>%
  group_by(Study, Brain_Region) %>%
  nest() %>%
  mutate(
    test = map(data, ~ cor.test(.x$values, .x$Power, method="spearman")), # S3 list-col
    tidied = map(test, tidy)
  ) %>%
  unnest(tidied) %>%
  select(-data, -test) %>%
  group_by(Study) %>%
  mutate(p_Bonferroni = p.adjust(p.value, method="bonferroni")) %>%
  arrange(desc(abs(estimate)))

# Plot regional SD correlations by brain region within the brain


# Bar plot of correlated regions
head_motion_SD_corr_regional %>%
  mutate(Brain_Region = fct_reorder(Brain_Region, abs(estimate))) %>%
  ggplot(data=., mapping=aes(x=Brain_Region,
                             y=abs(estimate),
                             fill=p_Bonferroni<0.05)) +
  geom_bar(stat="identity") +
  facet_grid(Study ~ ., scales="free", switch="both",
             space="free") +
  coord_flip() +
  scale_fill_manual(values=list("TRUE"="#4C7FC0",
                                "FALSE"="gray80")) +
  labs(fill = "Significant") +
  theme(strip.background = element_blank(),
        legend.position="bottom",
        strip.placement="outside")
ggsave(paste0(plot_path, "SD_Correlation_with_Movement_by_Region.svg"),
       width = 9, height = 16, units="in", dpi=300)

# Plot brain-wide average SD vs head movement
SD_brain_wide_avg <- UCLA_CNP_catch25 %>%
  plyr::rbind.fill(., ABIDE_ASD_catch25) %>%
  filter(names == "DN_Spread_Std") %>%
  group_by(Sample_ID) %>%
  summarise(meanval = mean(values, na.rm=T)) %>%
  left_join(., plyr::rbind.fill(UCLA_CNP_mean_FD, ABIDE_ASD_mean_FD)) %>%
  left_join(., plyr::rbind.fill(UCLA_CNP_sample_metadata, ABIDE_ASD_sample_metadata)) %>%
  mutate(Study = ifelse(Study == "UCLA_CNP", "UCLA CNP", "ABIDE"))  %>%
  mutate(Diagnosis = case_when(Diagnosis == "Schizophrenia" ~ "SCZ",
                               Diagnosis == "Bipolar" ~ "BP",
                               T ~ Diagnosis))

SD_brain_wide_avg %>%
  ggplot(data=., mapping=aes(x=Power, y=meanval, color=Diagnosis)) +
  geom_point(alpha=0.5) +
  facet_grid(Study ~ ., scale="free", switch="both") +
  ylab("Brain-wide BOLD SD") +
  xlab("Mean framewise displacement (FD)") +
  scale_color_manual(values = c("Control" = "#439E47", 
                                "SCZ" = "#9d60a8", 
                                "BP" = "#2F77C0", 
                                "ADHD" = "#e45075", 
                                "ASD" = "#E28328")) +
  stat_smooth(method="lm", se=F) +
  theme(strip.background = element_blank(),
        legend.position="bottom",
        strip.placement="outside") +
  stat_cor(method="spearman")
ggsave(paste0(plot_path, "mFD_vs_SD.svg"),
       width = 3.75, height=5.75, units="in", dpi=300)

# Calculate correlations
SD_brain_wide_avg %>%
  group_by(Study) %>%
  nest() %>%
  mutate(
    test = map(data, ~ cor.test(.x$meanval, .x$Power)), # S3 list-col
    tidied = map(test, tidy)
  ) %>%
  unnest(tidied)


################################################################################
# Compare balanced accuracy using just movement data for SVM
################################################################################

if (!file.exists(glue("{data_path}/UCLA_CNP_ABIDE_ASD_movement_SVM_results.feather"))) {
  
  # Import normalisation methods from python
  source_python("~/Library/CloudStorage/OneDrive-TheUniversityofSydney(Students)/github/fMRI_FeaturesDisorders/helper_functions/classification/mixed_sigmoid_normalisation.py")
  
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
  combined_movement_SVM_res <- do.call(plyr::rbind.fill, list(scz_mvmt_res, bpd_mvmt_res, adhd_mvmt_res, asd_mvmt_res)) %>%
    mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                        Comparison_Group == "Bipolar" ~ "BP",
                                        T ~ Comparison_Group)) %>%
    mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BP", "ADHD", "ASD")))
  
  # Save to feather
  pyarrow_feather$write_feather(combined_movement_SVM_res, glue("{data_path}/UCLA_CNP_ABIDE_ASD_movement_SVM_results.feather"))
} else {
  combined_movement_SVM_res <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_movement_SVM_results.feather"))
} 

combined_movement_SVM_res_summary <- combined_movement_SVM_res %>%
  group_by(Comparison_Group) %>%
  summarise(Balanced_Accuracy_Mean = mean(Balanced_Accuracy),
            Balanced_Accuracy_SD = sd(Balanced_Accuracy))

combined_movement_SVM_res %>%
  ggplot(data=., mapping=aes(x=Comparison_Group, y=100*Balanced_Accuracy)) +
  geom_violinhalf(aes(fill=Comparison_Group), scale="width",position = position_nudge(x=0.2))  +
  geom_boxplot(width=0.1, notch=FALSE, notchwidth = 0.4, outlier.shape = NA,
               fill=NA, color="white",
               position = position_nudge(x=0.27), coef = 0) +
  geom_hline(yintercept = 50, linetype=2, alpha=0.5) +
  geom_point(aes(color = Comparison_Group), position = position_jitterdodge(dodge.width = 1,
                                                                            jitter.width = 0.7),
             size = 1) +
  scale_fill_manual(values=c("Control" = "#439E47", 
                              "SCZ" = "#9d60a8", 
                              "BP" = "#2F77C0", 
                              "ADHD" = "#e45075", 
                              "ASD" = "#E28328")) +
  scale_color_manual(values=c("Control" = "#439E47", 
                             "SCZ" = "#9d60a8", 
                             "BP" = "#2F77C0", 
                             "ADHD" = "#e45075", 
                             "ASD" = "#E28328")) +
  ylab("Fold Balanced Accuracy (%)") +
  xlab("Clinical Group") +
  theme(legend.position = "none")
ggsave(paste0(plot_path, "SVM_movement_balacc.svg"),
       width = 3, height=3, units="in", dpi=300)

# Summarise across folds
combined_movement_SVM_res_summary %>%
  mutate(Balanced_Accuracy_Mean = 100*Balanced_Accuracy_Mean,
         Balanced_Accuracy_SD = 100*Balanced_Accuracy_SD)

combined_movement_SVM_res_summary


################################################################################

# p-values
# Functions to calculate empirical p-value
compare_main_and_null <- function(main_df_iter, null_distribution_df) {
  # Filter null to this data -- keep all grouping vars in this analysis type
  null_distribution_df <- null_distribution_df %>%
    dplyr::select(-index, -group_var, -Null_Iter_Number) %>%
    semi_join(., main_df_iter)
  
  # Compare main balanced accuracy with that of the empirical null distribution
  main_balanced_accuracy_value <- main_df_iter$Balanced_Accuracy_Mean
  null_balanced_accuracy_values <- null_distribution_df$Null_Balanced_Accuracy
  
  # Find proportion of iterations for which the main balanced accuracy is greater
  prop_main_greater <- sum(main_balanced_accuracy_value > null_balanced_accuracy_values)/length(null_balanced_accuracy_values)
  
  # Find p-value
  p_value <- 1 - prop_main_greater
  
  # Organise into dataframe
  main_df_iter$p_value <- p_value
  return(main_df_iter)
}

# Compare with the catch25 feature nulls
univariate_null_distribution <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_null_balanced_accuracy_distributions.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set,
         Analysis_Type == "Univariate_TS_Feature") %>%
  mutate(Comparison_Group = case_when(Comparison_Group == "Schizophrenia" ~ "SCZ",
                                      Comparison_Group == "Bipolar" ~ "BP",
                                      T ~ Comparison_Group))
         
combined_movement_SVM_res_summary_split <- split.data.frame(combined_movement_SVM_res_summary, combined_movement_SVM_res_summary$Comparison_Group)

movement_SVM_p_values <- combined_movement_SVM_res_summary_split %>%
  purrr::map_df(~ compare_main_and_null(main_df_iter = .x,
                                        null_distribution_df = univariate_null_distribution))