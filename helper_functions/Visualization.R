library(tidyverse)
library(patchwork)
library(scales)

#-------------------------------------------------------------------------------
# Plot distribution of classifier accuracies across given feature variable
# along with control vs schz proportions
#-------------------------------------------------------------------------------
plot_class_acc_w_props <- function(class_res,
                                   rdata_path,
                                   group_var = NULL,
                                   plot_title = "",
                                   noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"),
                                   ylab = "Number of ROIs") {
  
  # Calculate the proportion of control subjects after omitting NAs for each method
  ctrl_prop_list <- list()
  for (noise_proc in noise_procs) {
    # Clean up names
    noise_label <- gsub("\\+", "_", noise_proc)
    
    # Load corresponding feature matrix
    feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22_zscored.Rds", 
                                                         noise_label)))
    
    group_var_name <- ifelse(group_var == "Feature", "names", group_var)
    
    # Calculate proportion of controls after dropping NA
    if (!is.null(group_var)) {
      ctrl_proportion <- feature_matrix %>%
        group_by(Subject_ID, get(group_var_name)) %>%
        filter(!any(is.na(values))) %>%
        ungroup() %>%
        distinct(Subject_ID, group) %>%
        summarise(ctrl_prop = sum(group=="Control") / n()) %>%
        mutate(Noise_Proc = noise_proc)
    } else {
      ctrl_proportion <- feature_matrix %>%
        group_by(Subject_ID) %>%
        filter(!any(is.na(values))) %>%
        ungroup() %>%
        distinct(Subject_ID, group) %>%
        summarise(ctrl_prop = sum(group=="Control") / n()) %>%
        mutate(Noise_Proc = noise_proc)
    }
    ctrl_prop_list <- rlist::list.append(ctrl_prop_list, ctrl_proportion)
  }
  ctrl_prop <- do.call(plyr::rbind.fill, ctrl_prop_list)
  
  # Plot accuracy + balanced accuracy in histograms
  # Control subject proportion is highlighted for accuracy, 
  # 0.5 is highlighted for balanced accuracy
  if (!is.null(group_var)) {
    p_data <- class_res %>%
      dplyr::select(grouping_var, Sample_Type, Noise_Proc, accuracy, balanced_accuracy)
  } else {
    p_data <- class_res %>%
      dplyr::select(Noise_Proc, Sample_Type, accuracy, balanced_accuracy)
  }
  
  p_data %>%
    left_join(., ctrl_prop) %>%
    pivot_longer(cols=c(accuracy, balanced_accuracy),
                 names_to = "Metric",
                 values_to = "Value") %>%
    mutate(Full_Metric = ifelse(Sample_Type == "In-sample",
                                paste0("in_", Metric),
                                paste0("out_", Metric))) %>%
    mutate(Full_Metric = gsub("_", " ", Full_Metric),
           ctrl_prop = ifelse(Full_Metric %in% c("in balanced accuracy",
                                            "out balanced accuracy"), 
                              0.5, ctrl_prop)) %>%
    mutate(Full_Metric = factor(Full_Metric, levels = c("in accuracy", "out accuracy",
                                              "in balanced accuracy", 
                                              "out balanced accuracy")),
           Noise_Proc = factor(Noise_Proc, levels = noise_procs)) %>%
    ggplot(data=., mapping=aes(x=Value)) +
    geom_histogram(fill="lightsteelblue", bins=50) +
    ggtitle(plot_title) +
    geom_vline(aes(xintercept = ctrl_prop), linetype=2, color="gray30") +
    facet_grid(Noise_Proc ~ Full_Metric, scales="free_y", switch="y",
               labeller = labeller(Full_Metric = label_wrap_gen(20))) +
    xlab("Metric Value") +
    ylab(ylab) +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle=0),
          plot.title = element_text(hjust=0.5))
}


#-------------------------------------------------------------------------------
# Plot main vs null distribution histogram
#-------------------------------------------------------------------------------

plot_main_vs_null_hist <- function(main_res,
                                   null_res,
                                   xlab = "Value",
                                   ylab = "Scaled Density",
                                   title = "Main vs Null Distribution") {
  
  #extract hex color codes for a plot with three elements in ggplot2 
  hex_colors <- hue_pal()(2)

  null_res_long <- null_res %>%
    pivot_longer(cols = c(accuracy, balanced_accuracy),
                 names_to = "Metric",
                 values_to = "Values") %>%
    mutate(Metric = stringr::str_to_title(str_replace_all(Metric, "_", " ")),
           Noise_Proc = factor(Noise_Proc, levels = c("AROMA+2P",
                                                      "AROMA+2P+GMR",
                                                      "AROMA+2P+DiCER")))
  
  if (!("Sample_Type" %in% colnames(null_res_long))) {
    null_res_in <- null_res_long %>% mutate(Full_Metric = paste0("In ", Metric))
    null_res_out <- null_res_long %>% mutate(Full_Metric = paste0("Out ", Metric)) 
    null_res_long <- plyr::rbind.fill(null_res_in, null_res_out)%>%
      mutate(Full_Metric = factor(Full_Metric, levels = c("In Accuracy", "Out Accuracy",
                                                          "In Balanced Accuracy",
                                                          "Out Balanced Accuracy"))) 
  } else {
    null_res_long <- null_res_long %>%
      mutate(Full_Metric = ifelse(Sample_Type == "In-sample",
                                  paste0("In ", Metric),
                                  paste0("Out ", Metric)))%>%
      mutate(Full_Metric = factor(Full_Metric, levels = c("In Accuracy", "Out Accuracy",
                                                          "In Balanced Accuracy",
                                                          "Out Balanced Accuracy"))) 
  }
  
  main_res %>%
    dplyr::select(grouping_var, Sample_Type, Noise_Proc, accuracy, balanced_accuracy) %>%
    mutate(Type = "main") %>%
    pivot_longer(cols=c(accuracy, balanced_accuracy),
                 names_to = "Metric",
                 values_to = "Values") %>%
    mutate(Full_Metric = ifelse(Sample_Type == "In-sample",
                                paste0("in_", Metric),
                                paste0("out_", Metric))) %>%
    mutate(Full_Metric = stringr::str_to_title(str_replace_all(Full_Metric, "_", " "))) %>%
    mutate(Full_Metric = factor(Full_Metric, levels = c("In Accuracy", "Out Accuracy",
                                                        "In Balanced Accuracy",
                                                        "Out Balanced Accuracy"))) %>%
    mutate(Noise_Proc = factor(Noise_Proc, levels = c("AROMA+2P",
                                                      "AROMA+2P+GMR",
                                                      "AROMA+2P+DiCER"))) %>%
    ggplot(data=., mapping=aes(x=Values, fill = Type)) +
    scale_fill_manual(values = hex_colors) +
    geom_histogram(aes(y=0.5*..density..), 
                   bins = 50,
                   alpha=0.6, position="identity") +
    geom_histogram(data = null_res_long,
                   aes(y=0.5*..density..),
                   bins = 50,
                   alpha = 0.6,
                   position = "identity") +
    facet_grid(Noise_Proc ~ Full_Metric, switch="y", scales="fixed",
               labeller = labeller(Full_Metric = label_wrap_gen(20))) +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(title) +
    labs(fill = "Distribution") +
    theme(strip.text.y.left = element_text(angle=0),
          plot.title = element_text(hjust=0.5),
          strip.placement = "outside",
          legend.position = "bottom",
          legend.direction = "horizontal") 
}

#-------------------------------------------------------------------------------
# Pick top 6 grouping variables (e.g. brain region or feature) 
# in main vs null distributions, 
#-------------------------------------------------------------------------------

plot_top_6_vars_main_vs_null <- function(class_res_pvals,
                                         null_res,
                                         sample_type = "Out-of-sample",
                                         xlab = "Value",
                                         ylab = "Scaled Density",
                                         xloc = 0.6,
                                         yloc = 8.5,
                                         title = "Main vs Null Distribution",
                                         combo = FALSE) {
  
  # Filter null dataset if it contains a sample type column
  if ("Sample_Type" %in% colnames(null_res)) {
    null_res <- null_res %>% filter(Sample_Type == sample_type)
  }
  
  if (combo) {
    p <- class_res_pvals %>%
      ggplot(data=.) +
      geom_histogram(data = subset(null_res, Type=="null"),
                     aes(x = balanced_accuracy, y=0.5*..density..), 
                     bins = 50,
                     alpha=0.6, position="identity",
                     fill = "gray70") +
      geom_vline(mapping=aes(xintercept = balanced_accuracy, color = Noise_Proc), size=1.2) +
      ggtitle(title) +
      xlab(xlab) +
      ylab(ylab) +
      labs(color = "Noise Processing") +
      theme(strip.text.y.left = element_text(angle=0),
            strip.placement = "outside",
            legend.position = "bottom",
            legend.direction = "horizontal",
            plot.title = element_text(hjust=0.5))
  } else {
    top_features <- class_res_pvals %>%
      filter(Noise_Proc == "AROMA+2P",
             Sample_Type == sample_type) %>%
      arrange(desc(balanced_accuracy)) %>%
      slice(1:6) %>%
      pull(grouping_var)
    
    # Truncate p-value labels
    top_feature_plabs <- truncate_p_values(class_res_pvals, 3) %>%
      filter(Noise_Proc == "AROMA+2P", 
             grouping_var %in% top_features,
             Sample_Type == sample_type)
    
    p <- top_feature_plabs %>%
      mutate(grouping_var = factor(grouping_var, levels = top_features)) %>%
      ggplot(data=.) +
      geom_histogram(data = null_res %>% 
                       dplyr::select(balanced_accuracy, Noise_Proc) %>%
                       dplyr::filter(Noise_Proc == "AROMA+2P"),
                     aes(x=balanced_accuracy, y=0.5*..density..),
                     fill = "gray70", bins=50) +
      ggtitle(title) +
      geom_vline(aes(xintercept = balanced_accuracy), color = "red") +
      facet_wrap(grouping_var ~ ., scales="free_y", nrow = 2) +
      xlab(xlab) +
      ylab(ylab) +
      xlab("Balanced Accuracy") +
      geom_text(data = top_feature_plabs  %>%
                  mutate(grouping_var = factor(grouping_var, levels = top_features)),
                aes(label = paste0("P = ", bal_acc_p, "\nBH-FDR = ", bal_acc_p_adj)), 
                x = xloc, y = yloc) +
      theme(plot.title = element_text(hjust=0.5))
  }
  return(p)
  
}


#-------------------------------------------------------------------------------
# Visualise t-stat histograms for non-normalised and 
# normalised data across all ROIs
#-------------------------------------------------------------------------------

t_stat_histograms <- function(t_test_res, 
                              noise_proc = "AROMA+2P",
                              norm_methods = c("non-normalised",
                                               "z-score",
                                               "RobustSigmoid")) {
  noise_label <- gsub("\\+", "_", noise_proc)
  theme_set(cowplot::theme_cowplot())
  outer_plot_list <- list()
  
  t_stat_np <- subset(t_test_res, Noise_Proc==noise_proc)
  
  t_stat_limits <- c(min(t_stat_np$statistic, na.rm=T),
                     max(t_stat_np$statistic, na.rm=T))
  
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
        ggplot(data=., mapping=aes(x=statistic)) +
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
