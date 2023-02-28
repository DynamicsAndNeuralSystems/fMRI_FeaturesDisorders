################################################################################
# Load libraries and data
################################################################################
library(feather)
library(tidyverse)
library(glue)
library(icesTAF)
library(cowplot)
library(ggseg)
library(ggsegHO)
library(knitr)
library(kableExtra)
library(RColorBrewer)
library(correctR)
library(rstatix)
library(broom)
library(colorspace)
library(patchwork)
library(ggnewscale)
library(scales)
library(ggpubr)
library(ggsignif)
theme_set(theme_cowplot())

# Load in metadata 
load("data/schizophrenia_case_study/metadata.Rda")

# Load in movement data
load("data/schizophrenia_case_study/movement_data.Rda")

# Load in number of subjects retained per movement threhsold
load("data/schizophrenia_case_study/movement_threshold_counts.Rda")

# Load catch2 and catch22 feature values for data
load("data/schizophrenia_case_study/all_feature_values.Rda")

# Load balanced accuracy results from 10-resample 10-fold CV linear SVM
load("data/schizophrenia_case_study/CV_SVM_balanced_accuracy.Rda")

# Load SVM coefficients
load("data/schizophrenia_case_study/SVM_coefficients.Rda")

# Load null balanced accuracy results from 10-resample 10-fold CV linear SVM
load("data/schizophrenia_case_study/Null_SVM_balanced_accuracy.Rda")

# Load in univariate time-series feature info
TS_feature_info <- read.csv("info/catch24_info.csv")

# Load movement-based SVM results
load("data/schizophrenia_case_study/movement_based_SVM_results.Rda")

################################################################################
# Helper functions
################################################################################

# Functions to calculate empirical p-value
compare_main_and_null <- function(main_df_iter, null_distribution_df) {
  # Filter null to this data -- keep all grouping vars in this analysis type
  null_distribution_df <- null_distribution_df %>%
    dplyr::select(-index, -group_var) %>%
    semi_join(., main_df_iter)
  
  # Compare main balanced accuracy with that of the empirical null distribution
  main_balanced_accuracy_value <- main_df_iter$Balanced_Accuracy_Across_Repeats
  null_balanced_accuracy_values <- null_distribution_df$Null_Balanced_Accuracy
  
  # Find proportion of iterations for which the main balanced accuracy is greater
  prop_main_greater <- sum(main_balanced_accuracy_value > null_balanced_accuracy_values)/length(null_balanced_accuracy_values)
  
  # Find p-value
  p_value <- 1 - prop_main_greater
  
  # Organise into dataframe
  main_df_iter$p_value <- p_value
  return(main_df_iter)
}

################################################################################
# Compare catch2 vs. catch22 classification performance
################################################################################

# Aggregate the main results across all folds, independent of repeat
mean_balanced_accuracy <- SVM_balanced_accuracy_all_folds %>%
  group_by(Univariate_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
  
  # First get the average across folds within a repeat
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T)) %>%
  
  # Then get the average across repeats
  group_by(Univariate_Feature_Set, Analysis_Type, group_var) %>%
  summarise(Balanced_Accuracy_Across_Repeats = mean(Balanced_Accuracy_Across_Folds))

# Summary for SD
SVM_balanced_accuracy_all_folds %>%
  filter(group_var == "Combo") %>%
  group_by(Repeat_Number, Univariate_Feature_Set) %>%
  summarise(mean_balacc = mean(Balanced_Accuracy)) %>%
  ungroup() %>%
  group_by(Univariate_Feature_Set) %>%
  summarise(meanb = mean(mean_balacc),
            sdb = sd(mean_balacc))

# Split the dataframe and calculate empirical p-values
balanced_accuracy_split <- mean_balanced_accuracy %>%
  group_by(Analysis_Type, Univariate_Feature_Set, group_var) %>%
  group_split()

p_values <- balanced_accuracy_split %>%
  purrr::map_df(~ compare_main_and_null(main_df_iter = .x,
                                        null_distribution_df = null_balanced_accuracy_all_folds)) %>%
  group_by(Analysis_Type, Univariate_Feature_Set) %>%
  mutate(p_value_BH = p.adjust(p_value, method="BH"),
         p_value_Bonferroni = p.adjust(p_value, method="bonferroni"))


# Brain region-wise: plot regions that are significantly better/worse with catch2 relative to catch22
mean_balanced_accuracy %>%
  filter(Analysis_Type == "Univariate_Brain_Region") %>%
  left_join(., p_values) %>%
  group_by(Analysis_Type, group_var) %>%
  mutate(Status = case_when(Balanced_Accuracy_Across_Repeats[Univariate_Feature_Set=="catch2"] - Balanced_Accuracy_Across_Repeats[Univariate_Feature_Set=="catch22"] > 0.05 ~  "catch2",
                            Balanced_Accuracy_Across_Repeats[Univariate_Feature_Set=="catch22"] - Balanced_Accuracy_Across_Repeats[Univariate_Feature_Set=="catch2"] > 0.05 ~  "catch22",
                            T ~ "Neither")) %>%
  filter(Status != "Neither") %>%
  mutate(important_p_value = p_value_Bonferroni[Univariate_Feature_Set==Status]) %>%
  filter(important_p_value < 0.05) %>%
  mutate(Status = paste0(Status, " better")) %>%
  ggplot(data = ., mapping=aes(x = Univariate_Feature_Set, 
                               y=Balanced_Accuracy_Across_Repeats, 
                               group=group_var, 
                               color = Status)) +
  scale_colour_brewer(palette = "Dark2") +
  ggtitle("Schizophrenia vs. Control Classification\nwithin Brain Regions by Feature Set") +
  ylab("Mean Balanced Accuracy") +
  xlab("Feature Set") +
  facet_grid(. ~ Status) +
  geom_point() +
  geom_line(alpha=0.8) +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5))
ggsave(glue("output/Schizophrenia_BalAcc_Brain_Region_Comparison_5percent.png"),
       width = 6, height = 3.5, units="in", dpi=300)

# TS Feature-wise
mean_balanced_accuracy %>%
  filter(Analysis_Type == "Univariate_TS_Feature") %>%
  left_join(., p_values) %>%
  filter(p_value_Bonferroni < 0.05) %>%
  ggplot(data=., mapping=aes(x = Balanced_Accuracy_Across_Repeats, color = Univariate_Feature_Set)) +
  geom_vline(aes(xintercept = Balanced_Accuracy_Across_Repeats, color = Univariate_Feature_Set), linewidth=1.5) +
  labs(color = "Feature Set") +
  scale_colour_brewer(palette = "Dark2") +
  ggtitle("Schizophrenia vs. Control Classification within TS Features by Feature Set") +
  xlab("Mean Balanced Accuracy") +
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        plot.title = element_text(hjust=0.5,margin=margin(0,0,20,0),
                                  size=11))
ggsave("output/Schizophrenia_BalAcc_TSFeature_Feature_Set.png",
       width = 6, height = 2, units="in", dpi=300)

# Combo-wise
SVM_balanced_accuracy_all_folds %>%
  group_by(Univariate_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
  
  # First get the average across folds within a repeat
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T)) %>%
  
  filter(Analysis_Type == "Univariate_Combo") %>%
  ggplot(data=., mapping=aes(x = Univariate_Feature_Set, 
                             y = Balanced_Accuracy_Across_Folds,
                             fill = Univariate_Feature_Set)) +
  ggtitle("Schizophrenia vs. Control\nClassification with\nRegion-by-Feature\nCombinations") +
  scale_fill_brewer(palette = "Dark2") +
  labs(fill = "Feature Set") +
  geom_violin() +
  geom_line(aes(group = Repeat_Number), alpha=0.3) +
  geom_boxplot(fill=NA, width=0.1) +
  ylab("Mean Balanced Accuracy by Resample") +
  xlab("Feature Set") +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5, size=13))
ggsave("output/Schizophrenia_BalAcc_Combo_Feature_Set.png",
       width = 4, height = 6, units="in", dpi=300)

# Use correctR to test for difference across repeated k-fold CV
num_samples <- length(unique(metadata$Sample_ID))

SVM_balanced_accuracy_all_folds %>%
  filter(Analysis_Type == "Univariate_Combo") %>%
  dplyr::rename("model" = "Univariate_Feature_Set",
                k = "Fold",
                r = "Repeat_Number",
                values = "Balanced_Accuracy") %>%
  repkfold_ttest(data = .,
                 n1 = ceiling(0.9*num_samples),
                 n2 = floor(0.1*num_samples),
                 k = 10,
                 r = 10)

################################################################################
# Visualize standard deviation + mean in the brain
################################################################################

# Helper function to run t-test for given statistic
run_t_test_for_feature <- function(all_feature_values, input_feature) {
  results <- all_feature_values %>%
    filter(names==input_feature) %>%
    dplyr::select(Brain_Region, Diagnosis, values) %>%
    mutate(Diagnosis = factor(Diagnosis, levels = c("Schizophrenia", "Control"))) %>%
    group_by(Brain_Region) %>%
    nest() %>%
    mutate(
      fit = map(data, ~ t.test(values ~ Diagnosis, data = .x)),
      tidied = map(fit, tidy)
    ) %>% 
    unnest(tidied) %>%
    dplyr::select(-data, -fit) %>%
    arrange(p.value) %>%
    ungroup() %>%
    mutate(label = ifelse(str_detect(Brain_Region, "ctx-"),
                          gsub("-", "_", Brain_Region),
                          as.character(Brain_Region))) %>%
    mutate(label = gsub("ctx_", "", label)) %>%
    dplyr::select(label, statistic)
  return(results)
}

extract_coefs_for_feature <- function(SVM_coefficients, input_feature) {
  feature_coefs <- SVM_coefficients %>%
    filter(Analysis_Type=="Univariate_TS_Feature",
           group_var==input_feature) %>%
    arrange(desc(abs(coef))) %>%
    mutate(label = ifelse(str_detect(`Feature Name`, "ctx-"),
                          gsub("-", "_", `Feature Name`),
                          as.character(`Feature Name`))) %>%
    mutate(label = gsub("ctx_", "", label)) %>%
    dplyr::select(label, coef) %>%
    dplyr::rename("statistic" = "coef")
  
  return(feature_coefs)
}

# Run t-test for SD
SD_Tdata_for_ggseg <- run_t_test_for_feature(all_feature_values = all_feature_values,
                                             input_feature = "DN_Spread_Std")

# Run t-test for Mean
Mean_Tdata_for_ggseg <- run_t_test_for_feature(all_feature_values = all_feature_values,
                                               input_feature = "DN_Mean")

# Helper function to plot cortex data
plot_data_in_cortex <- function(results_data, min_value, 
                                max_value, 
                                label_title) {
  cortex_plot <- results_data %>%
    ggseg(atlas = dk, mapping = aes(fill = statistic),
          position = "stacked", colour = "gray40") +
    scale_fill_continuous_divergingx(palette = 'RdBu', 
                                     mid = 0, 
                                     rev=TRUE, 
                                     na.value="gray90",
                                     limits = c(min_value,
                                                max_value)) +
    theme_void() +
    labs(fill=label_title) +
    guides(fill = guide_colorbar(title.position = "top", 
                                 nrow = 1,
                                 barwidth = 10, 
                                 barheight = 0.75,
                                 title.hjust = 0.5,
                                 label.position = "bottom"))  +
    theme(plot.title = element_blank(),
          legend.position = "bottom") 
  
  return(cortex_plot)
}

plot_data_in_subcortex <- function(results_data, min_value, max_value, label_title) {
  subcortex_plot <- results_data %>%
    ggplot() +
    geom_brain(atlas = aseg, mapping = aes(fill=statistic), 
               side = "coronal", colour = "gray40") +
    scale_fill_continuous_divergingx(palette = 'RdBu', 
                                     mid = 0, 
                                     rev=TRUE, 
                                     na.value="gray90",
                                     limits = c(min_value,
                                                max_value)) +
    theme_void() +
    labs(fill=label_title) +
    guides(fill = guide_colorbar(title.position = "top", 
                                 nrow = 1,
                                 barwidth = 10, 
                                 barheight = 0.75,
                                 title.hjust = 0.5,
                                 label.position = "bottom"))  +
    theme(plot.title = element_blank(),
          legend.position = "bottom") 
  
  
  return(subcortex_plot)
}

# Plot SD T-statistics in the brain
SD_T_cortex_plot <- plot_data_in_cortex(results_data = SD_Tdata_for_ggseg,
                                        min_value = min(SD_Tdata_for_ggseg$statistic),
                                        max_value = max(SD_Tdata_for_ggseg$statistic),
                                        label_title = "Region T-statistic for SD")

SD_T_subcortex_plot <- plot_data_in_subcortex(results_data = SD_Tdata_for_ggseg,
                                              min_value = min(SD_Tdata_for_ggseg$statistic),
                                              max_value = max(SD_Tdata_for_ggseg$statistic),
                                              label_title = "Region T-statistic for SD")

# Wrap the two plots and combine legends
wrap_plots(list(SD_T_cortex_plot, SD_T_subcortex_plot),
           nrow = 2, heights = c(0.65,0.35)) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom',
        legend.title = element_text(size=9),
        legend.justification = "center")
ggsave("output/Schizophrenia_SD_Tstats_Cortex_Subcortex.png",
       width = 3, height = 4, units="in", dpi=300)


# Plot SD SVM coefficients in the brain
SD_coefs_for_ggseg <- extract_coefs_for_feature(SVM_coefficients, input_feature="DN_Spread_Std")

SD_coefs_cortex_plot <- plot_data_in_cortex(results_data = SD_coefs_for_ggseg,
                                            min_value = min(SD_coefs_for_ggseg$statistic),
                                            max_value = max(SD_coefs_for_ggseg$statistic),
                                            label_title = "Region SVM Coefficients for SD")

SD_coefs_subcortex_plot <- plot_data_in_subcortex(results_data = SD_coefs_for_ggseg,
                                                  min_value = min(SD_coefs_for_ggseg$statistic),
                                                  max_value = max(SD_coefs_for_ggseg$statistic),
                                                  label_title = "Region SVM Coefficients for SD")
# Wrap the two plots and combine legends
wrap_plots(list(SD_coefs_cortex_plot, SD_coefs_subcortex_plot),
           nrow = 2, heights = c(0.65,0.35)) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom',
        legend.title = element_text(size=9),
        legend.justification = "center")
ggsave("output/Schizophrenia_SD_Coefs_Cortex_Subcortex.png",
       width = 3, height = 4, units="in", dpi=300)


# Plot t-statistic for Mean in schizophrenia vs controls by region
Mean_T_cortex_plot <- plot_data_in_cortex(results_data = Mean_Tdata_for_ggseg,
                                          min_value = min(Mean_Tdata_for_ggseg$statistic),
                                          max_value = max(Mean_Tdata_for_ggseg$statistic),
                                          label_title = "Region T-statistic for Mean")

Mean_T_subcortex_plot <- plot_data_in_subcortex(results_data = Mean_Tdata_for_ggseg,
                                                min_value = min(Mean_Tdata_for_ggseg$statistic),
                                                max_value = max(Mean_Tdata_for_ggseg$statistic),
                                                label_title = "Region T-statistic for Mean")
# Wrap the two plots and combine legends
wrap_plots(list(Mean_T_cortex_plot, Mean_T_subcortex_plot),
           nrow = 1, widths = c(0.7,0.3)) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom',
        legend.justification = "center")
ggsave("output/Schizophrenia_Mean_Tstats_Cortex_Subcortex.png",
       width = 4, height = 3, units="in", dpi=300)


# SVM Coefficients for Mean
Mean_coefs_for_ggseg <- extract_coefs_for_feature(SVM_coefficients, input_feature="DN_Mean")

Mean_coefs_cortex_plot <- plot_data_in_cortex(results_data = Mean_coefs_for_ggseg,
                                              min_value = min(Mean_coefs_for_ggseg$statistic),
                                              max_value = max(Mean_coefs_for_ggseg$statistic),
                                              label_title = "Region SVM Coefficients for Mean")

Mean_coefs_subcortex_plot <- plot_data_in_subcortex(results_data = Mean_coefs_for_ggseg,
                                                    min_value = min(Mean_coefs_for_ggseg$statistic),
                                                    max_value = max(Mean_coefs_for_ggseg$statistic),
                                                    label_title = "Region SVM Coefficients for Mean")
# Wrap the two plots and combine legends
wrap_plots(list(Mean_coefs_cortex_plot, Mean_coefs_subcortex_plot),
           nrow = 1, widths = c(0.7,0.3)) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom',
        legend.justification = "center")
ggsave("output/Schizophrenia_Mean_Coefs_Cortex_Subcortex.png",
       width = 4, height = 3, units="in", dpi=300)


# Wrap the two plots and combine legends
wrap_plots(list(Mean_coefs_cortex_plot, Mean_coefs_subcortex_plot),
           nrow = 1, widths = c(0.7,0.3)) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom',
        legend.justification = "center")
ggsave("output/Schizophrenia_Mean_Coefs_Cortex_Subcortex.png",
       width = 4, height = 3, units="in", dpi=300)

################################################################################
# Movement analysis
################################################################################

# Plot mFD Power by Diagnosis
mFD_data %>%
  ggplot(data=., mapping=aes(x=Diagnosis, y=mFD_Power, fill=Diagnosis)) +
  scale_fill_brewer(palette = "Dark2") +
  geom_violin() +
  geom_boxplot(fill=NA, width=0.1) +
  ggtitle("Mean Framewise Displacement\nby Diagnosis") +
  ylab("mFD-Power") +
  geom_hline(yintercept = 0.55, alpha=0.7) +
  geom_hline(yintercept = 0.25, alpha=0.7) +
  geom_signif(test = "wilcox.test",
              comparisons = list(c("Schizophrenia", "Control")), 
              map_signif_level=TRUE) +
  theme(legend.position="none",
        plot.title = element_text(hjust=0.5, size=12))
ggsave("output/Schizophrenia_Movement_Distributions.png", width=4, height=4, 
       units="in", dpi=300)

# Summarise mFD in control vs schizophrenia
mFD_data %>%
  group_by(Diagnosis) %>%
  summarise(mean_mFD = round(mean(mFD_Power), 2),
            SD_mFD = round(sd(mFD_Power), 2))

# Is there a significant difference?
mFD_data %>%
  wilcox_test(mFD_Power ~ Diagnosis)

# Plot number of subjects retained at each threshold
movement_threshold_counts <- movement_threshold_counts %>%
  mutate(Threshold = factor(Threshold, levels=c("None", 
                                                "Lenient", 
                                                "Stringent"))) 

movement_threshold_counts %>%
  ggplot(data=., 
         mapping=aes(x=Diagnosis, y=Threshold)) +
  geom_tile(color="white", data=subset(movement_threshold_counts, 
                                       Diagnosis=="Control"),
            aes(fill=n)) +
  scale_fill_gradientn(colors = c(alpha(brewer.pal(4, "Dark2")[1], 0.7),
                                  brewer.pal(4, "Dark2")[1])) +
  new_scale_fill() +
  geom_tile(color="white", data=subset(movement_threshold_counts, 
                                       Diagnosis=="Schizophrenia"),
            aes(fill=n)) +
  scale_fill_gradientn(colors = c(alpha(brewer.pal(4, "Dark2")[2], 0.7),
                                  brewer.pal(4, "Dark2")[2])) +
  geom_text(aes(label=n)) +
  ggtitle("Number of Participants Retained\nper Head Movement Threshold") +
  ylab("Threshold Type") +
  scale_y_discrete(limits=rev) +
  theme(legend.position="none",
        plot.title = element_text(hjust=0.5, size=12))
ggsave("output/Schizophrenia_Movement_Num_Subjects.png", width=4, height=2.5, 
       units="in", dpi=300)

# Compare performance of features at no, lenient, and stringent movement thresholds
averaged_movement_SVM_res <- all_movement_SVM_res %>%
  group_by(Threshold_Type, Univariate_Feature_Set, Analysis_Type, group_var, Repeat_Number) %>%
  summarise(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy)) %>%
  group_by(Threshold_Type, Univariate_Feature_Set, Analysis_Type, group_var) %>%
  summarise(Balanced_Accuracy_Across_Repeats = mean(Balanced_Accuracy_Across_Folds))

# Line plots
averaged_movement_SVM_res %>%
  left_join(., p_values %>% dplyr::select(-Balanced_Accuracy_Across_Repeats)) %>%
  filter(p_value_Bonferroni < 0.05) %>%
  mutate(Threshold_Type = str_replace_all(Threshold_Type, " Movement Threshold", ""),
         Analysis_Type = case_when(Analysis_Type == "Univariate_Brain_Region" ~ "Region-wise",
                                   Analysis_Type == "Univariate_TS_Feature" ~ "Feature-wise",
                                   Analysis_Type == "Univariate_Combo" ~ "Combo-wise")) %>%
  mutate(Threshold_Type = ifelse(Threshold_Type=="No", "None", Threshold_Type)) %>%
  mutate(Threshold_Type = factor(Threshold_Type, levels=c("None", 
                                                          "Lenient", 
                                                          "Stringent")),
         Analysis_Type = factor(Analysis_Type, levels=c("Region-wise", "Feature-wise", "Combo-wise"))) %>%
  mutate(plot_label = paste(Analysis_Type, Univariate_Feature_Set)) %>%
  ggplot(data=., mapping=aes(x=Threshold_Type,
                             y=Balanced_Accuracy_Across_Repeats, 
                             group = group_var,
                             color = plot_label)) +
  geom_line(alpha=0.6, linewidth=0.8) +
  scale_colour_brewer(palette = "Dark2") +
  ggtitle("TS Feature-wise Balanced Accuracy\nby Movement Threshold") +
  ylab("Balanced Accuracy Across Repeats") +
  xlab("Movement Threshold Type") +
  facet_grid(Analysis_Type ~ Univariate_Feature_Set, scales="free",
             switch="y", space="free") +
  theme(legend.position="none",
        plot.title = element_text(hjust=0.5, size=13),
        axis.text.x = element_text(size=10),
        strip.placement = "outside")
ggsave(glue("output/Schizophrenia_Movement_SVM_results.png"),
       width = 5, height = 6, units="in", dpi=300, bg="white")

# Compare mFD values with feature values
feature_motion_lm <- all_feature_values %>%
  mutate(Age = as.numeric(Age)) %>%
  left_join(., mFD_data) %>%
  group_by(names, Brain_Region) %>%
  nest() %>%
  mutate(
    test = map(data, ~ lm(values ~ mFD_Power + Diagnosis + Age, data=.)), # S3 list-col
    tidied = map(test, tidy)
  ) %>% 
  unnest(tidied) %>%
  dplyr::select(-data, -test) %>%
  ungroup() %>%
  filter(term != "(Intercept)") %>%
  group_by(term) %>%
  mutate(p_adj = p.adjust(p.value, method="bonferroni")) %>%
  filter(p_adj < 0.05)

# Plot regions for which mFD_Power is still significantly associated with SD
mFD_power_lm_plot <- feature_motion_lm %>%
  filter(term == "mFD_Power") %>%
  mutate(names = case_when(names == "DN_Spread_Std" ~ "SD",
                           T ~ str_replace_all(names, "_", " "))) %>%
  group_by(names) %>%
  mutate(Brain_Region = fct_reorder(Brain_Region, estimate, .fun=max)) %>%
  ungroup() %>%
  mutate(names = fct_reorder(names, estimate, .fun=mean, .desc=T)) %>%
  ggplot(data=., mapping=aes(x = Brain_Region, y=estimate)) +
  geom_bar(stat="identity", aes(fill=estimate)) +
  scale_fill_gradientn(colors = c(alpha(brewer.pal(4, "Dark2")[1], 0.5),
                                  brewer.pal(4, "Dark2")[1])) +
  facet_grid(names ~ ., scales="free", space="free", switch="y",
             labeller = labeller(names = label_wrap_gen(20))) +
  xlab("Brain Region") +
  ylab("Coefficient Estimate") +
  coord_flip() +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle=0),
        legend.position = "none")

age_lm_plot <- feature_motion_lm %>%
  filter(term == "Age") %>%
  mutate(label = formatC(estimate, format = "e", digits=1)) %>%
  mutate(label = ifelse(abs(estimate) < 0.00001, label, NA),
         names = case_when(names == "DN_Spread_Std" ~ "SD",
                           names == "DN_Mean" ~ "Mean",
                           T ~ str_replace_all(names, "_", " "))) %>%
  group_by(names) %>%
  mutate(Brain_Region = fct_reorder(Brain_Region, estimate, .fun=max)) %>%
  ungroup() %>%
  mutate(names = fct_reorder(names, estimate, .fun=mean, .desc=T)) %>%
  ggplot(data=., mapping=aes(x = Brain_Region, y=estimate)) +
  geom_bar(stat="identity", aes(fill=estimate)) +
  scale_fill_gradientn(colors = c(brewer.pal(4, "Dark2")[2]),
                       alpha(brewer.pal(4, "Dark2")[2], 0.5)) +
  xlab("Brain Region") +
  ylab("Coefficient Estimate") +
  facet_grid(names ~ ., scales="free", space="free", switch="y",
             labeller = labeller(names = label_wrap_gen(15))) +
  coord_flip() +
  geom_text(aes(label = label),
            color = brewer.pal(4, "Dark2")[2],
            hjust=1) +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle=0),
        legend.position = "none")

wrap_plots(list(mFD_power_lm_plot,
                age_lm_plot),
           heights=c(0.7, 0.3),
           ncol=1) + 
  plot_annotation(title = "Linear Regression Coefficients\nby Brain Region") &
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y = element_text(size=10))
ggsave("output/Schizophrenia_Movement_Age_LM_Coefficients.png",
       width = 6, height = 8, units="in", dpi=300)

# Look at SD
feature_motion_lm %>%
  filter(names == "DN_Spread_Std", term=="mFD_Power", p_adj<0.05) %>%
  arrange(desc(estimate))

