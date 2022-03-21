#------------------------------------
# This script contains functions for visualizations used in the readme
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 20 March 2022
#--------------------------------------

library(tidyverse)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Plot classification accuracy or t-statistic across all regions
plot_ROI_acc_by_feature <- function(region_wise_univ_class_res, 
                                    xlab="T statistic",
                                    noise_proc,
                                    plot_path) {
  theme_set(cowplot::theme_cowplot())
  # Clean up labels
  xlab_file <- gsub(" ", "_", xlab)
  noise_label <- gsub("\\+", "_", noise_proc)
  
  region_wise_univ_class_res %>%
    # filter(feature %in% unique(region_wise_univ_class_res$feature)[1:2]) %>%
    mutate(feature = str_replace_all(feature, "_", " ")) %>%
    mutate(feature = str_replace_all(feature, "catch22 ", "")) %>%
    ggplot(data=., mapping=aes(x=statistic_value)) +
    geom_histogram(fill="lightsteelblue") +
    geom_vline(xintercept=0, linetype=2) +
    xlab(xlab) +
    ylab("Number of ROIs") +
    facet_wrap(feature ~ ., scales="free", nrow=4,
               labeller = labeller(feature = label_wrap_gen(26)))
  ggsave(paste0(plot_path, sprintf("UCLA_%s_catch22_%s_histograms.png", 
                                   noise_label, xlab)),
         width=15, height=10, units="in", dpi=300)
}

