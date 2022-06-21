library(tidyverse)
library(patchwork)
library(scales)

#-------------------------------------------------------------------------------
# Plot distribution of classifier accuracies across given feature variable
# along with control vs schz proportions
#-------------------------------------------------------------------------------
plot_class_acc_w_props <- function(class_res,
                                   feature_set = "catch22",
                                   rdata_path,
                                   group_var = NULL,
                                   plot_title = "",
                                   noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"),
                                   ylab = "Number of ROIs") {
  
  # Calculate the proportion of control subjects after omitting NAs
  ctrl_prop <- readRDS(paste0(rdata_path, "Filtered_subject_info_",
                                 feature_set, ".Rds")) %>%
    dplyr::summarise(ctrl_prop = sum(group=="Control") / n()) %>%
    pull(ctrl_prop)
  
  # Plot accuracy + balanced accuracy in histograms
  # Control subject proportion is highlighted for accuracy, 
  # 0.5 is highlighted for balanced accuracy
  if (!is.null(group_var)) {
    p_data <- class_res %>%
      dplyr::select(grouping_var, Sample_Type, Noise_Proc, accuracy, balanced_accuracy)
  } else {
    p_data <- class_res %>%
      dplyr::select(Noise_Proc, Sample_Type, accuracy, accuracy_SD, balanced_accuracy, balanced_accuracy_SD)
  }
  
  p_data <- p_data %>%
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
           Noise_Proc = factor(Noise_Proc, levels = noise_procs)) 
  
  # Return histogram if the group var is provided
  if (!is.null(group_var)) {
    p_data %>%
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
  } else {
    p_data %>%
      pivot_longer(cols = c(accuracy_SD, balanced_accuracy_SD),
                   names_to = "Metric2", values_to = "SD") %>%
      mutate(Metric2 = str_replace_all(Metric2, "_SD", "")) %>%
      filter(Metric == Metric2) %>%
      dplyr::select(-Metric2) %>%
      ggplot(data=., mapping=aes(x=Full_Metric, y=Value)) +
      geom_bar(aes(fill = Full_Metric), stat = 'identity') +
      geom_errorbar(aes(ymin = Value - SD, ymax = Value + SD)) +
      ggtitle(plot_title) +
      facet_grid(Noise_Proc ~ ., scales="free_y", switch="y") +
      xlab("Performance Metric") +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
      ylab(ylab) +
      theme(strip.placement = "outside",
            strip.text.y.left = element_text(angle=0),
            plot.title = element_text(hjust=0.5),
            legend.position = 'none')
  }
  
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
    mutate(Metric = stringr::str_to_title(str_replace_all(Metric, "_", " ")))
  
  if ("Noise_Proc" %in% colnames(null_res_long)) {
    null_res_long <- null_res_long %>%
      dplyr::mutate(Noise_Proc = factor(Noise_Proc, levels = c("AROMA+2P",
                                                               "AROMA+2P+GMR",
                                                               "AROMA+2P+DiCER")))
  }
  
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
      filter(Sample_Type == sample_type) %>%
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
    
    # Prep null data
    if ("Noise_Proc" %in% colnames(null_res)) {
      null_data_for_plot <- null_res %>%
        dplyr::select(balanced_accuracy, Noise_Proc) %>%
        dplyr::filter(Noise_Proc == "AROMA+2P")
    } else {
      null_data_for_plot <- null_res %>%
        dplyr::select(balanced_accuracy)
    }
      
    p <- top_feature_plabs %>%
      mutate(grouping_var = factor(grouping_var, levels = top_features)) %>%
      ggplot(data=.) +
      geom_histogram(data = null_data_for_plot,
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
# Plot density distribution of main vs null results given a dataset
#-------------------------------------------------------------------------------

plot_main_vs_null_bal_acc_density <- function(main_res, null_res, pvals,
                                              grouping_type, result_color,
                                              line_only = FALSE) {

  
  # Density plot
  p <- main_res %>%
    filter(Sample_Type == "Out-of-sample",
           Noise_Proc == "AROMA+2P+GMR") %>%
    left_join(., pvals) %>%
    ungroup() %>%
    ggplot(data=., mapping=aes(x = balanced_accuracy)) +
    ggtitle(grouping_type) +
    xlab("10-Fold CV Balanced Accuracy") +
    ylab("Density") 
  
  if (line_only) {
    raw_bal_acc_vector <- null_res %>%
      dplyr::filter(Sample_Type == "Out-of-sample" &
                      Noise_Proc == "AROMA+2P+GMR") %>%
      pull(balanced_accuracy)
    
    raw_p_thresh <- quantile(raw_bal_acc_vector, probs = c(0.95))
    
    p <- p +
      geom_density(data = subset(null_res,
                                 Sample_Type == "Out-of-sample" &
                                   Noise_Proc == "AROMA+2P+GMR"),
                   aes(fill = "Null"),
                   alpha = 0.7) +
      geom_vline(aes(xintercept = balanced_accuracy,
                     fill = "Main"),
                 color = result_color, size = 1.6) +
      geom_vline(aes(xintercept = raw_p_thresh),
                 color = "black", size = 1.2,
                 linetype = 2)
  } else {
    # Find cutoff value for BH-adjusted significance for group-wise results
    group_wise_bal_acc_threshold <- main_res %>%
      filter(Sample_Type == "Out-of-sample",
             Noise_Proc == "AROMA+2P+GMR") %>%
      left_join(., pvals) %>%
      ungroup() %>%
      arrange(bal_acc_p_adj) %>%
      mutate(is_sig = bal_acc_p_adj < 0.05) %>%
      distinct(is_sig, .keep_all = T) %>%
      filter(!is_sig)
    
    p <- p +
      geom_density(aes(fill = "Main")) +
      geom_density(data = subset(null_res,
                                 Sample_Type == "Out-of-sample" &
                                   Noise_Proc == "AROMA+2P+GMR"),
                   aes(fill = "Null"),
                   alpha = 0.7) +
      geom_vline(data = group_wise_bal_acc_threshold,
                 mapping = aes(xintercept = balanced_accuracy),
                 linetype = 2, size=1.2)
  }

  p <- p +
    labs(fill = "Result") +
    scale_fill_manual(values = c(result_color, "gray40")) +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5))
  
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
