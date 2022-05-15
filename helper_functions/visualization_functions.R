library(tidyverse)

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
