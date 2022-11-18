################################################################################
# Load libraries
################################################################################

library(tidyverse)
library(icesTAF)
library(cowplot)
library(theft)
library(ggridges)
theme_set(theme_cowplot())

################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/FigureS5/")
icesTAF::mkdir(plot_path)

ASD_data_path <- "~/data/ABIDE_ASD/"
ASD_rdata_path <- paste0(ASD_data_path, "processed_data/Rdata/")

ROI_univariate_nulls <- readRDS(paste0(ASD_rdata_path,
                                       "ABIDE_ASD_ROI_wise_model_permutation_null_catch22_inv_prob.Rds"))

ROI_univariate_nulls %>%
  ggplot(data=., mapping=aes(x=balanced_accuracy, fill=grouping_var)) +
  geom_histogram() +
  ggtitle("ASD -- Raw") +
  facet_wrap(grouping_var ~ ., scales="free_y", ncol=7,
             labeller = labeller(grouping_var = label_wrap_gen(28))) +
  xlab("Balanced Accuracy") +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5),
        strip.text = element_text(size=10),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave(paste0(plot_path, "ASD_catch22_ROI_wise_null_dist.png"),
       width = 15, height = 9, units="in", dpi=300)

ROI_univariate_nulls %>%
  ggplot(data=., mapping=aes(x=balanced_accuracy, y=grouping_var,
                             fill=grouping_var)) +
  geom_density_ridges2(stat = "binline", alpha=0.8) +
  theme(legend.position="none") +
  ylab("Harvard-Oxford Cortical Brain Region") +
  xlab("Null Balanced Accuracy")
ggsave(paste0(plot_path, "ASD_catch22_ROI_wise_null_dist.png"),
       width = 9, height = 10, units="in", dpi=300)