#------------------------------------
# This script contains general visualization functions
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 20 March 2022
#--------------------------------------

library(tidyverse)
library(cowplot)
library(janitor)
library(patchwork)
library(ComplexHeatmap)
library(ggpubr)
library(patchwork)
library(circlize)
theme_set(theme_cowplot())

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Plot feature-derived PC1 vs. PC2 for a given ROI

# Generic dimplot function
dimplot_helper <- function(data, data_wide, title, subtitle) {
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
  p <- fits %>%
    ggplot(aes(x = .fitted1, y = .fitted2)) +
    geom_point(size = 1.5, aes(colour = group)) +
    ggplot2::labs(title = title,
                  subtitle = subtitle,
                  x = paste0("PC 1"," (",eigen_pc1,")"),
                  y = paste0("PC 2"," (",eigen_pc2,")"),
                  colour = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "bottom",
                   legend.text = element_text(size=12),
                   plot.title = element_text(hjust=0.5)) + 
    guides(colour = guide_legend(override.aes = list(size=3)))
  
  return(p)
}

# Main PCA dimplot function
PCA_dimplot_for_ROI <- function(this_ROI, rdata_path, 
                                noise_procs =  c("AROMA+2P",
                                                 "AROMA+2P+GMR",
                                                 "AROMA+2P+DiCER"),
                                norm_methods = c("z-score",
                                                 "RobustSigmoid")) {
  
  # Clean brain region label
  region_label <- this_ROI %>%
    gsub("ctx-lh-", "Left ", .) %>%
    gsub("ctx-rh-", "Right ", .) 
  
  # Instantiate list for plots
  plot_list <- list()
  
  # Iterate over each noise processing method
  for (noise_proc in noise_procs) {
    
    noise_label <- gsub("\\+", "_", noise_proc)
    
    # Non-normalized data
    feature_matrix <- readRDS(paste0(rdata_path, "UCLA_", 
                                     noise_label, "_catch22.Rds"))
    
    
    feature_lookup <- data.frame(names = unique(feature_matrix$names)) %>%
      mutate(feature_clean = janitor::make_clean_names(names))
    
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
    
    title = noise_proc
    subtitle = "Non-normalised data"
    
    # Call helper dimplot function
    noise_proc_plot <- dimplot_helper(data = ROI_mat,
                                      data_wide = ROI_data_wide,
                                      title = title,
                                      subtitle = subtitle)
    
    plot_list <- rlist::list.append(plot_list, noise_proc_plot)
    
    ############################################################################
    for (norm_method in norm_methods) {

      # Run feature normalisation
      normed <- normalise_feature_frame(feature_matrix, 
                                        names_var = "names", 
                                        values_var = "values", 
                                        method = norm_method)
      
      # Filter the catch22 dataset to the given ROI
      ROI_data <- normed %>%
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
      
      title = noise_proc
      subtitle = paste0(norm_method, " normalised")
      
      # Call helper dimplot function
      noise_proc_plot <- dimplot_helper(data = ROI_mat,
                                        data_wide = ROI_data_wide,
                                        title = title,
                                        subtitle = subtitle)
      
      plot_list <- rlist::list.append(plot_list, noise_proc_plot)
    }
  }
  
  p <- patchwork::wrap_plots(plot_list, ncol = 3) + 
    plot_annotation(
      title = region_label,
      theme = theme(plot.title = element_text(hjust=0.5))
    ) + 
    plot_layout(guides = "collect") & 
    theme(legend.position = 'bottom')
  
  return(p)
}


#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Visualise non-normalised and normalised data for a given ROI with violin plots
violin_plots_for_ROI <- function(this_ROI, 
                                 rdata_path, 
                                 noise_proc = "AROMA+2P",
                                 norm_methods = c("z-score",
                                                  "RobustSigmoid")) {
  noise_label <- gsub("\\+", "_", noise_proc)
  
  # Clean brain region label
  region_label <- this_ROI %>%
    gsub("ctx-lh-", "Left ", .) %>%
    gsub("ctx-rh-", "Right ", .) 
  
  # Iterate over each noise processing method
  ROI_data_list <- list()
  
  # Non-normalized data
  feature_matrix <- readRDS(paste0(rdata_path, "UCLA_", 
                                   noise_label, "_catch22.Rds")) %>%
    mutate(group = factor(group, levels=c("Control", "Schz")))
  
  
  feature_lookup <- data.frame(names = unique(feature_matrix$names)) %>%
    mutate(feature_clean = janitor::make_clean_names(names))
  
  # Filter the catch22 dataset to the given ROI
  ROI_data <- feature_matrix %>%
    filter(Brain_Region==this_ROI) %>%
    mutate(Norm_Method = "non-normalised",
           Noise_Proc = noise_proc)
  
  # Append data to list
  ROI_data_list <- rlist::list.append(ROI_data_list, ROI_data)
  
  for (norm_method in norm_methods) {
    
    feature_matrix_norm <- normalise_feature_frame(feature_matrix, 
                                                   names_var = "names", 
                                                   values_var = "values", 
                                                   method = norm_method)
    
    # Filter the catch22 dataset to the given ROI
    ROI_data <- feature_matrix_norm %>%
      filter(Brain_Region==this_ROI) %>%
      mutate(Norm_Method = norm_method,
             Noise_Proc = noise_proc)
    
    # Append data to list
    ROI_data_list <- rlist::list.append(ROI_data_list, ROI_data)
  }
  
  # Combine data
  ROI_data_full <- do.call(plyr::rbind.fill, ROI_data_list) 
  
  outer_plot_list <- list()
  for (group in c(1,2)) {
    
    if (group==1) {
      ROI_data_group <- ROI_data_full %>%
        filter(names %in% unique(ROI_data_full$names)[1:11])  
    } else {
      ROI_data_group <- ROI_data_full %>%
        filter(names %in% unique(ROI_data_full$names)[12:22])  
    }
    
    ROI_t_test <- ROI_data_group %>%
      group_by(names, Noise_Proc, Norm_Method) %>%
      nest() %>%
      mutate(
        test = map(data, ~ t.test(.x$values ~ .x$group)), # S3 list-col
        tidied = map(test, tidy)
      ) %>% 
      unnest(tidied) %>%
      dplyr::select(-data, -test) %>%
      dplyr::select(names, Noise_Proc, Norm_Method, estimate) %>%
      mutate(estimate = round(estimate, 3)) %>%
      dplyr::rename("t_stat" = "estimate") %>%
      left_join(., ROI_data_full) %>%
      group_by(names, Norm_Method) %>%
      mutate(max_value = max(values, na.rm=T)) %>%
      distinct(names, Noise_Proc, Norm_Method, t_stat, max_value)
    
    plot_list <- list()
    for (norm_method in c("non-normalised", norm_methods)) {
      data_subset <- subset(ROI_data_group, Norm_Method==norm_method)
      
      # Create base plot
      p <- data_subset %>%
        left_join(., ROI_t_test) %>%
        mutate(names = str_replace_all(names, "_", " ")) %>%
        ggplot(data=., mapping=aes(x=group, y=values)) +
        geom_violin(aes(fill=group)) +
        geom_boxplot(width=0.25, fill=NA) +
        geom_text(aes(x=1.5, y=max_value*1.1, label=t_stat),
                  size=5) +
        ggtitle(norm_method) +
        scale_y_continuous(expand=c(0,0,0.1,0)) +
        facet_grid(names ~ ., scales = "free", switch="y",
                   labeller = labeller(names = label_wrap_gen(20)))
      
      if (norm_method=="non-normalised") {
        p <- p +
          theme(strip.text.y.left = element_text(angle=0),
                strip.placement = "outside",
                axis.text.x = element_blank(),
                axis.title.x=element_blank(),
                axis.ticks.x = element_blank(),
                legend.position="bottom",
                legend.direction="horizontal",
                plot.title = element_text(hjust=0.5))
        
      } else {
        p <- p +
          theme(strip.text.y.left = element_blank(),
                strip.background = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_blank(),
                axis.title.x=element_blank(),
                axis.ticks.x = element_blank(),
                legend.position="bottom",
                legend.direction="horizontal",
                plot.title = element_text(hjust=0.5))
        
      }
      
      # Append to list
      plot_list <- rlist::list.append(plot_list, p)
    }
    
    p <- patchwork::wrap_plots(plot_list, 
                          ncol = length(plot_list)) + 
      plot_layout(guides = "collect") & 
      theme(legend.position = 'bottom')
    outer_plot_list <- rlist::list.append(outer_plot_list, p)
  }
  final_p <- wrap_plots(outer_plot_list, ncol = length(outer_plot_list)) + 
    plot_annotation(title = sprintf("Distributions for %s %s catch22 features",
                                    region_label, noise_proc),
                    theme = theme(plot.title = element_text(size = 20))) + 
    plot_layout(guides = "collect") & 
    theme(legend.position = 'bottom')
  print(final_p)
}


#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Visualise t-stat histograms for non-normalised and 
# normalised data across all ROIs
t_stat_histograms <- function(t_test_res, 
                              noise_proc = "AROMA+2P",
                              norm_methods = c("non-normalised",
                                               "z-score",
                                               "RobustSigmoid")) {
  noise_label <- gsub("\\+", "_", noise_proc)
  theme_set(cowplot::theme_cowplot())
  outer_plot_list <- list()
  
  t_stat_np <- subset(t_test_res, Noise_Proc==noise_proc)
  
  t_stat_limits <- c(min(t_stat_np$statistic_value, na.rm=T),
                     max(t_stat_np$statistic_value, na.rm=T))
  
  for (group in c(1,2)) {
    
    if (group==1) {
      t_stat_group <- t_stat_np %>%
        filter(feature %in% unique(t_stat_np$feature)[1:11])  
    } else {
      t_stat_group <- t_stat_np %>%
        filter(feature %in% unique(t_stat_np$feature)[12:22])  
    }
    
    # Iterate over each normalisation method
    plot_list <- list()
    
    for (norm_method in norm_methods) {
      t_stat_subset <- subset(t_stat_group, Norm_Method==norm_method)
      
      # Create base plot
      p <- t_stat_subset %>%
        mutate(feature = str_replace_all(feature, "_", " ")) %>%
        ggplot(data=., mapping=aes(x=statistic_value)) +
        geom_histogram(fill="lightsteelblue") +
        geom_vline(xintercept=0, linetype=2) +
        scale_x_continuous(limits = t_stat_limits) +
        ggtitle(norm_method) +
        ylab("Number of ROIs") +
        xlab("T statistic") +
        facet_grid(feature ~ ., scales="free", switch="y",
                   labeller = labeller(feature = label_wrap_gen(20)))
    
      # Add or hide facet strips accordingly
      if (norm_method=="non-normalised") {
        p <- p  +
          theme(strip.placement = "outside",
                strip.text.y.left = element_text(angle=0),
                plot.title=element_text(hjust=0.5))
      } else {
        p <- p +
          theme(strip.placement = "outside",
                strip.text.y.left = element_blank(),
                strip.background = element_blank(),
                axis.title.y = element_blank(),
                plot.title=element_text(hjust=0.5))
      }
      
      plot_list <- rlist::list.append(plot_list, p)
    }
    
    p <- patchwork::wrap_plots(plot_list, 
                               ncol = length(plot_list)) 
    
    for (i in c(1, 3:length(plot_list))) {
      # Remove title from second subplot
      p[[i]] = p[[i]] + 
        theme(axis.title.x = element_blank() )
    }
    
    
    outer_plot_list <- rlist::list.append(outer_plot_list, p)
  }
  final_p <- wrap_plots(outer_plot_list, ncol = length(outer_plot_list))
  print(final_p)
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


#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Plot a hierarchically-clustered heatmap of feature X region showing 
# Univariate classification statistic
plot_class_stat_heatmap <- function(classification_results,
                                    statistic,
                                    noise_proc,
                                    norm_method,
                                    min_val=-4,
                                    max_val=4) {
  classification_results_wide <- classification_results %>%
    mutate(Brain_Region = gsub("ctx-lh-", "Left ", Brain_Region)) %>%
    mutate(Brain_Region = gsub("ctx-rh-", "Right ", Brain_Region)) %>%
    pivot_wider(id_cols=c("feature"),
                names_from = "Brain_Region",
                values_from = "statistic_value")
  classification_mat <- classification_results_wide %>%
    mutate(feature = gsub("catch22_", "", feature)) %>%
    column_to_rownames("feature") %>%
    as.matrix()
  
  col_func = circlize::colorRamp2(seq(min_val, max_val, length = 3), c("blue", "#EEEEEE", "red"))
  
  htmp <- ComplexHeatmap::Heatmap(classification_mat, 
                                  col = col_func,
                                  column_title = sprintf("Classification metrics for UCLA %s %s Data",
                                                         noise_proc, norm_method),
                                  heatmap_legend_param = list(title=statistic,
                                                              title_position = "leftcenter",
                                                              legend_width = unit(5, "cm"),
                                                              legend_direction="horizontal"),
                                  row_names_gp = gpar(fontsize = 10),
                                  column_names_gp = gpar(fontsize = 10))
  
  draw(htmp, heatmap_legend_side="bottom")
}

