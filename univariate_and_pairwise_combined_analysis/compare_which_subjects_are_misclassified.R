# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--project_path", default="/project/hctsa/annie/")
parser$add_argument("--github_dir", default="/project/hctsa/annie/github/")
parser$add_argument("--data_path", default="/project/hctsa/annie/data/UCLA_Schizophrenia/")
parser$add_argument("--rdata_path", default="/project/hctsa/annie/data/UCLA_Schizophrenia/Rdata/")
parser$add_argument("--feature_set", default="catch22")
# project_path <- "D:/Virtual_Machines/Shared_Folder/github/"
# github_dir <- "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/"
# data_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/UCLA_Schizophrenia/"
# rdata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/UCLA_Schizophrenia/Rdata/"
# feature_set <- "catch22"

# Parse input arguments
args <- parser$parse_args()
project_path <- args$project_path
github_dir <- args$github_dir
data_path <- args$data_path
rdata_path <- args$rdata_path
feature_set <- args$feature_set

### Source functions
# Main
source(paste0(github_dir, "helper_functions/Linear_SVM.R"))
source(paste0(github_dir, "helper_functions/Visualization.R"))
source(paste0(github_dir, "helper_functions/Null_distributions.R"))

set.seed(127)

# Compare all three noise processing methods
noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")

# Use e1071 SVM with a linear kernel
test_package = "e1071"
kernel = "linear"

# Load subject metadata
subject_metadata <- read.csv(paste0(data_path, "participants.csv")) %>%
  dplyr::rename("Subject_ID" = "sampleID")

# Univariate ROI-wise
univariate_roi_p_vals <- readRDS(paste0(rdata_path, "ROI_wise_CV_linear_SVM_model_permutation_null_catch22_inv_prob_pvals.Rds"))
univariate_roi_sig <- univariate_roi_p_vals %>%
  filter(bal_acc_p_adj < 0.05, 
         Noise_Proc == "AROMA+2P+GMR") %>%
  pull(grouping_var)

univariate_roi_subject_class <- readRDS(paste0(rdata_path, "ROI_wise_CV_linear_SVM_catch22_inv_prob.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR",
         grouping_var %in% univariate_roi_sig) %>%
  mutate(grouping_var = str_replace_all(grouping_var, "ctx-rh-|Right-", "Right "),
         grouping_var = str_replace_all(grouping_var, "ctx-lh-|Left-", "Left ")) %>%
  group_by(Subject_ID) %>%
  mutate(num_incorr = sum(!Prediction_Correct)) %>%
  ungroup() 

# Univariate feature-wise
univariate_feature_p_vals <- readRDS(paste0(rdata_path, "Feature_wise_CV_linear_SVM_model_permutation_null_catch22_inv_prob_pvals.Rds"))
univariate_feature_sig <- univariate_feature_p_vals %>%
  filter(bal_acc_p_adj < 0.05, 
         Noise_Proc == "AROMA+2P+GMR") %>%
  pull(grouping_var)

univariate_feature_subject_class <- readRDS(paste0(rdata_path, "Feature_wise_CV_linear_SVM_catch22_inv_prob.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR",
         grouping_var %in% univariate_feature_sig) %>%
  group_by(Subject_ID) %>%
  mutate(num_incorr = sum(!Prediction_Correct)) %>%
  ungroup() 

# Univariate combo-wise
univariate_combo_subject_class <- readRDS(paste0(rdata_path, "Combo_wise_CV_linear_SVM_catch22_inv_prob.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR") %>%
  mutate(grouping_var = "Combo") %>%
  group_by(Subject_ID) %>%
  mutate(num_incorr = sum(!Prediction_Correct)) %>%
  ungroup() 

# Merge data
merged_data <- do.call(plyr::rbind.fill,
                       list(univariate_roi_subject_class,
                            univariate_feature_subject_class,
                            univariate_combo_subject_class))

# Find # times each subject is misclassified across all univariate models
misclassifications_w_info <- merged_data %>%
  group_by(Subject_ID, Actual_Diagnosis) %>%
  summarise(num_incorr = sum(!Prediction_Correct)) %>%
  arrange(desc(num_incorr)) %>%
  left_join(., subject_metadata)

# Print the subjects with >=5 misclassifications
misclassifications_w_info %>%
  filter(num_incorr >= 5) %>%
  dplyr::select(Subject_ID:gender) %>%
  dplyr::select(-diagnosis) %>%
  knitr::kable() %>%
  kableExtra::kable_styling(full_width = F)

# By age
misclassifications_w_info %>%
  ggplot(data=., mapping=aes(x=age, y=num_incorr, color=diagnosis)) +
  geom_point() +
  ylab("# Incorrect Predictions") +
  xlab("Age") +
  theme(legend.position = "bottom")

# Violin plot of prediction correctness by group
misclassifications_w_info %>%
  ggplot(data=., mapping=aes(x=Actual_Diagnosis, y=num_incorr)) +
  geom_violin(aes(fill=Actual_Diagnosis)) +
  geom_boxplot(fill=NA, color="black", width=0.1) +
  xlab("Diagnosis") +
  ylab("# Incorrect Predictions") +
  theme(legend.position="bottom")

# Bar chart of correctness by group per significant feature
merged_data %>%
  group_by(Actual_Diagnosis, grouping_var) %>%
  summarise(correct_prop = sum(Prediction_Correct) / n()) %>%
  ggplot(data=., mapping=aes(y=grouping_var, x=correct_prop, fill=Actual_Diagnosis)) +
  geom_bar(stat="identity", position='dodge') +
  ylab("Linear SVM Feature") +
  xlab("% Correct") +
  scale_fill_manual(values = c("green", "red")) +
  labs(fill="") +
  theme(legend.position = "bottom")

# Heatmap of prediction correctness
guidebar <- merged_data %>%
  mutate(Subject_ID = fct_reorder(Subject_ID, num_incorr, .fun=max)) %>%
  arrange(Subject_ID) %>%
  ggplot(data = ., mapping=aes(x=Subject_ID, y=0, fill=Actual_Diagnosis)) +
  geom_tile() +
  theme_void() +
  scale_fill_manual(values = c("green", "red")) +
  theme(legend.position = "none")

correct_plot <- merged_data %>%
  mutate(Subject_ID = fct_reorder(Subject_ID, num_incorr, .fun=max)) %>%
  arrange(Subject_ID) %>%
  ggplot(data=., mapping=aes(x = Subject_ID, 
                             y = grouping_var,
                             fill = Prediction_Correct)) +
  geom_tile() +
  ylab("Linear SVM Feature") +
  labs(fill = "Correctly Predicted") +
  scale_fill_manual(values = c("black", "white")) +
  xlab("Subjects") +
  theme(axis.text.x = element_blank(),
        legend.position = "bottom")

guidebar / correct_plot + 
  plot_layout(heights = c(1, 40))
 