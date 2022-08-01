# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--project_path", default="/project/hctsa/annie/")
parser$add_argument("--github_dir", default="/project/hctsa/annie/github/")
parser$add_argument("--rdata_path", default="/project/hctsa/annie/data/scz/UCLA/Rdata/")
parser$add_argument("--feature_set", default="catch22")
# project_path <- "D:/Virtual_Machines/Shared_Folder/github/"
# github_dir <- "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/"
# rdata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/Rdata/"
# feature_set <- "catch22"

# Parse input arguments
args <- parser$parse_args()
project_path <- args$project_path
github_dir <- args$github_dir
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

# Violin plot of prediiction correctness by group
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
 