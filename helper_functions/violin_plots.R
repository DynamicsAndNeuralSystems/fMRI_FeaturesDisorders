#------------------------------------
# Functions to generate violin plots showing distribution of values by feature and/or ROI
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 20 March 2022
#--------------------------------------

library(tidyverse)
library(cowplot)
library(janitor)
theme_set(theme_cowplot())

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Plot distribution of schizophrenia vs. control values for a given ROI for the top N features
violin_plot_for_ROI <- function(feature_matrix, class_res, 
                                    this_ROI, num_feature=8) {
  
  feature_lookup <- data.frame(names = unique(feature_matrix$names)) %>%
    mutate(feature_clean = janitor::make_clean_names(names))
  
  # Clean brain region label
  region_label <- this_ROI %>%
    gsub("ctx-lh-", "Left ", .) %>%
    gsub("ctx-rh-", "Right ", .) 
  
  # Filter the catch22 dataset to the given ROI
  ROI_data <- feature_matrix %>%
    filter(Brain_Region==this_ROI)
  
  # Filter the univariate classification results to the given feature
  ROI_tstats <- UCLA_AROMA_2P_catch22_ROIwise_t_test %>%
    filter(Brain_Region == this_ROI)
  
  # Find the top n features
  top_features <- ROI_tstats %>%
    mutate(abs_t = abs(statistic_value)) %>%
    arrange(desc(abs_t)) %>%
    top_n(num_feature, abs_t) %>%
    pull(feature)
  
  # Make violin plot
  ROI_data %>%
    left_join(., feature_lookup) %>%
    mutate(feature = paste(method, feature_clean, sep="_")) %>%
    filter(feature %in% top_features) %>%
    mutate(feature = gsub("catch22_", "", feature)) %>%
    ggplot(data=., mapping=aes(x=group, y=values)) +
    geom_violin(aes(fill=group), alpha=0.5, color="gray20") +
    geom_boxplot(fill=NA, width=0.25, outlier.shape = NA) +
    facet_wrap(feature ~ ., scale="free_y") +
    ggtitle(sprintf("Top %s Features for %s",
                    num_feature, region_label)) +
    ylab("Raw Value") +
    labs(fill = "Diagnosis") +
    theme(legend.position="bottom",
          plot.title = element_text(hjust=0.5),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())
}


#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Plot distribution of schizophrenia vs. control values for a given feature for the top N ROIs
violin_plot_for_feature <- function(feature_matrix, class_res, 
                                    this_feature, num_ROI=8) {
  
  this_feature_t <- paste0("catch22_", make_clean_names(this_feature))
  
  # Filter the catch22 dataset to the given feature
  feature_data <- feature_matrix %>%
    filter(names==this_feature)
  
  # Filter the univariate classification results to the given feature
  feature_tstats <- UCLA_AROMA_2P_catch22_ROIwise_t_test %>%
    filter(feature == this_feature_t)
  
  # Find the top n regions
  top_regions <- feature_tstats %>%
    mutate(abs_t = abs(statistic_value)) %>%
    arrange(desc(abs_t)) %>%
    top_n(num_ROI, abs_t) %>%
    pull(Brain_Region)
  
  # Make violin plot
  feature_data %>%
    filter(Brain_Region %in% top_regions) %>%
    mutate(Brain_Region = gsub("ctx-lh-", "Left ", Brain_Region)) %>%
    mutate(Brain_Region = gsub("ctx-rh-", "Right ", Brain_Region)) %>%
    ggplot(data=., mapping=aes(x=group, y=values)) +
    geom_violin(aes(fill=group), alpha=0.5, color="gray20") +
    geom_boxplot(fill=NA, width=0.25, outlier.shape = NA) +
    facet_wrap(Brain_Region ~ ., scale="free_y") +
    ggtitle(sprintf("Top %s ROIs for %s",
                    num_ROI, this_feature_t)) +
    ylab("Raw Value") +
    labs(fill = "Diagnosis") +
    theme(legend.position="bottom",
          plot.title = element_text(hjust=0.5),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())
}