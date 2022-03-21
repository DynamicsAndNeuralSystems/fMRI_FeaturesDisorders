#------------------------------------
# This function plots PC1 vs. PC2 for a given feature and/or ROI
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
# Generic dimplot function
dimplot_helper <- function(data, data_wide, subtitle) {
  fits <- data %>%
    stats::prcomp(scale = FALSE)
  
  # Retrieve eigenvalues and tidy up variance explained for plotting
  eigens <- fits %>%
    broom::tidy(matrix = "eigenvalues") %>%
    dplyr::filter(PC %in% c(1,2)) %>% # Filter to just the 2 going in the plot
    dplyr::select(c(PC, percent)) %>%
    dplyr::mutate(percent = round(percent*100), digits = 1)
  
  eigen_pc1 <- eigens %>%
    dplyr::filter(PC == 1)
  
  eigen_pc2 <- eigens %>%
    dplyr::filter(PC == 2)
  
  eigen_pc1 <- paste0(eigen_pc1$percent,"%")
  eigen_pc2 <- paste0(eigen_pc2$percent,"%")
  
  fits <- fits %>%
    broom::augment(data) %>%
    dplyr::rename(id = `.rownames`) %>%
    dplyr::mutate(id = as.factor(id)) %>%
    dplyr::rename(.fitted1 = .fittedPC1,
                  .fitted2 = .fittedPC2) %>%
    mutate(group = data_wide$group)
  
  # Make PCA dimplot
  fits %>%
    ggplot(aes(x = .fitted1, y = .fitted2)) +
    geom_point(size = 1.5, aes(colour = group)) +
    ggplot2::labs(title = "Low dimensional projection of time-series",
                  subtitle = subtitle,
                  x = paste0("PC 1"," (",eigen_pc1,")"),
                  y = paste0("PC 2"," (",eigen_pc2,")"),
                  colour = NULL) +
    ggplot2::scale_colour_brewer(palette = "Dark2") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "bottom")
}

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Plot feature-derived PC1 vs. PC2 for a given ROI
PCA_dimplot_for_ROI <- function(feature_matrix, class_res, 
                                    this_ROI) {
  
  feature_lookup <- data.frame(names = unique(feature_matrix$names)) %>%
    mutate(feature_clean = janitor::make_clean_names(names))
  
  # Clean brain region label
  region_label <- this_ROI %>%
    gsub("ctx-lh-", "Left ", .) %>%
    gsub("ctx-rh-", "Right ", .) 
  
  # Filter the catch22 dataset to the given ROI
  ROI_data <- feature_matrix %>%
    filter(Brain_Region==this_ROI)
  
  # Run PCA
  ROI_data_wide <- ROI_data %>%
    left_join(., feature_lookup) %>%
    mutate(feature = paste(method, feature_clean, sep="_")) %>%
    group_by(feature, group) %>%
    mutate(values = ifelse(is.na(values), median(values, na.rm=T), values)) %>%
    pivot_wider(id_cols=c(Subject_ID, group), 
                names_from=feature,
                values_from=values)
  
  ROI_mat <- ROI_data_wide %>%
    ungroup() %>%
    select(-Subject_ID, -group) %>%
    as.matrix()
  
  subtitle = sprintf("Brain region: %s", region_label)
  
  # Call helper dimplot function
  dimplot_helper(data = ROI_mat,
                 data_wide = ROI_data_wide,
                 subtitle = subtitle)
  
}

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Plot region-derived PC1 vs. PC2 for a given feature
PCA_dimplot_for_feature <- function(feature_matrix, class_res, 
                                this_feature) {
  
  feature_lookup <- data.frame(names = unique(feature_matrix$names)) %>%
    mutate(feature_clean = janitor::make_clean_names(names))
  
  this_feature_label <- subset(feature_lookup, names==this_feature) %>% 
    pull(feature_clean)
  
  # Filter the catch22 dataset to the given ROI
  feature_data <- feature_matrix %>%
    filter(names==this_feature)
  
  # Run PCA
  feature_data_wide <- feature_data %>%
    left_join(., feature_lookup) %>%
    mutate(feature = paste(method, feature_clean, sep="_")) %>%
    group_by(feature, group) %>%
    mutate(values = ifelse(is.na(values), median(values, na.rm=T), values)) %>%
    pivot_wider(id_cols=c(Subject_ID, group), 
                names_from=Brain_Region,
                values_from=values)
  
  feature_mat <- feature_data_wide %>%
    ungroup() %>%
    select(-Subject_ID, -group) %>%
    as.matrix()
  
  subtitle = sprintf("Feature: %s", this_feature_label)
  
  # Call helper dimplot function
  dimplot_helper(data = feature_mat,
                 data_wide = feature_data_wide,
                 subtitle = subtitle)
  
}
