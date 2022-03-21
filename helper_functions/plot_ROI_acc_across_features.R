#------------------------------------
# This script contains functions for visualizations used in the readme
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 20 March 2022
#--------------------------------------

library(tidyverse)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Plot classification accuracy or t-statistic across features for a given region
plot_feature_acc_by_ROI <- function(class_res, 
                                    xlab="T statistic",
                                    noise_proc,
                                    plot_path) {
  theme_set(cowplot::theme_cowplot())
  # Clean up labels
  xlab_file <- gsub(" ", "_", xlab)
  noise_label <- gsub("\\+", "_", noise_proc)
  
  class_res %>%
    mutate(Brain_Region = gsub("ctx-lh-", "Left ", Brain_Region)) %>%
    mutate(Brain_Region = gsub("ctx-rh-", "Right ", Brain_Region)) %>%
    ggplot(data=., mapping=aes(x=statistic_value)) +
    geom_histogram(fill="lightsteelblue") +
    geom_vline(xintercept=0, linetype=2) +
    xlab(xlab) +
    ylab("Number of Features") +
    facet_wrap(Brain_Region ~ ., scales="free", ncol=13,
               labeller = labeller(Brain_Region = label_wrap_gen(26)))
  ggsave(paste0(plot_path, sprintf("UCLA_%s_catch22_%s_feature_wise_histograms.png", 
                                   noise_label, xlab)),
         width=22, height=10, units="in", dpi=300)
}

