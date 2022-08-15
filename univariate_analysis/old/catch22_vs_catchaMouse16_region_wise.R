project_path <- "D:/Virtual_Machines/Shared_Folder/github/"
github_dir <- "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/"
rdata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/Rdata/"

### Source functions
# Main
source(paste0(github_dir, "helper_functions/Linear_SVM.R"))
source(paste0(github_dir, "helper_functions/Visualization.R"))
source(paste0(github_dir, "helper_functions/Null_distributions.R"))
library(patchwork)

set.seed(127)

# Use AROMA+2P+GMR noise-processing
noise_proc = "AROMA+2P+GMR"
noise_label = "AROMA_2P_GMR"


# Load data for catch22 and catchaMouse16
catch22_data <- readRDS(paste0(rdata_path, "ROI_wise_CV_linear_SVM_catch22_inv_prob.Rds")) %>%
  filter(Noise_Proc == noise_proc,
         Sample_Type == "Out-of-sample")
catchaMouse16_data <- readRDS(paste0(rdata_path, "ROI_wise_CV_linear_SVM_catchaMouse16_inv_prob.Rds")) %>%
  filter(Noise_Proc == noise_proc,
         Sample_Type == "Out-of-sample")

# Delta bar chart
delta_plot <- plyr::rbind.fill(catch22_data, catchaMouse16_data) %>%
  pivot_wider(id_cols = c(grouping_var),
              names_from = feature_set,
              values_from = balanced_accuracy) %>%
  rowwise() %>%
  mutate(delta_catchaMouse16_catch22 = catchaMouse16 - catch22,
         max_bacc = max(c(catch22, catchaMouse16))) %>%
  ungroup() %>%
  mutate(grouping_var = str_replace_all(grouping_var, "ctx-lh-|Left-", "Left "),
         grouping_var = str_replace_all(grouping_var, "ctx-rh-|Right-", "Right ")) %>%
  mutate(grouping_var = fct_reorder(grouping_var, max_bacc, .fun=max)) %>%
  arrange(grouping_var) %>%
  ggplot(data=., mapping=aes(x = delta_catchaMouse16_catch22,
                             y = grouping_var,
                             fill = delta_catchaMouse16_catch22)) +
  geom_bar(stat="identity") +
  geom_vline(xintercept = 0, linetype = 2) + 
  xlab("\u0394 Balanced Accuracy with catchaMouse16 vs. catch22") +
  ylab("Brain Region") +
  theme(legend.position = "none") +
  scale_fill_gradient2(low = "#2166AC", mid = "gray90", high = "#B2182B")

# Heatmap
bacc_plot <- plyr::rbind.fill(catch22_data, catchaMouse16_data) %>%
  mutate(grouping_var = str_replace_all(grouping_var, "ctx-lh-|Left-", "Left "),
         grouping_var = str_replace_all(grouping_var, "ctx-rh-|Right-", "Right ")) %>%
  mutate(grouping_var = fct_reorder(grouping_var, balanced_accuracy, .fun=max)) %>%
  ggplot(data=., mapping=aes(x = feature_set, 
                             y = grouping_var,
                             fill = balanced_accuracy)) +
  geom_tile() +
  geom_text(aes(label = round(100*balanced_accuracy, 1)),
            color = "black")  +
  scale_fill_gradient(low = "white", high = "deepskyblue2") +
  xlab("Feature Set") +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

delta_plot + bacc_plot