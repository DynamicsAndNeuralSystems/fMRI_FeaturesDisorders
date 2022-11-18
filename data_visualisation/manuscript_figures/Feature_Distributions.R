################################################################################
# Load libraries
################################################################################

library(tidyverse)
library(icesTAF)
library(cowplot)
library(theft)
theme_set(theme_cowplot())

################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
source("Manuscript_Draft_Visualisations_Helper.R")
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/FigureS2/")
icesTAF::mkdir(plot_path)

ASD_data_path <- "~/data/ABIDE_ASD/"
ASD_rdata_path <- paste0(ASD_data_path, "processed_data/Rdata/")

################################################################################
# Figure A: univariate raw distributions
################################################################################

# ASD raw data
ASD_univariate_catch22_raw <- readRDS(paste0(ASD_rdata_path, "ABIDE_ASD_catch22_filtered.Rds"))
ASD_univariate_catch22_raw %>%
  mutate(names = str_replace_all(names, "_", " ")) %>%
  ggplot(data=., mapping = aes(x=values, fill=names)) +
  geom_histogram() +
  ggtitle("ASD -- Raw") +
  facet_wrap(names ~ ., scales="free", ncol=5,
             labeller = labeller(names = label_wrap_gen(22))) +
  xlab("Feature Values") +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave(paste0(plot_path, "ASD_catch22_raw.png"),
       width = 10, height = 9, units="in", dpi=300)

# ASD z-scored data
ASD_univariate_catch22_zscore <- readRDS(paste0(ASD_rdata_path, "ABIDE_ASD_catch22_filtered_zscored.Rds"))
ASD_univariate_catch22_zscore %>%
  mutate(names = str_replace_all(names, "_", " ")) %>%
  ggplot(data=., mapping = aes(x=values, fill=names)) +
  geom_histogram() +
  ggtitle("ASD -- z-scored") +
  facet_wrap(names ~ ., scales="free", ncol=5,
             labeller = labeller(names = label_wrap_gen(25))) +
  xlab("Feature Values") +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave(paste0(plot_path, "ASD_catch22_zscored.png"),
       width = 10, height = 9, units="in", dpi=300)

# ASD robust sigmoid transformed data
ASD_univariate_catch22_robustsigmoid <- normalise_feature_frame(data = ASD_univariate_catch22_raw,
                                                                method = "RobustSigmoid")
ASD_univariate_catch22_robustsigmoid %>%
  mutate(names = str_replace_all(names, "_", " ")) %>%
  ggplot(data=., mapping = aes(x=values, fill=names)) +
  geom_histogram() +
  ggtitle("ASD -- Robust Sigmoid") +
  facet_wrap(names ~ ., scales="free_y", ncol=5,
             labeller = labeller(names = label_wrap_gen(25))) +
  xlab("Feature Values") +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave(paste0(plot_path, "ASD_catch22_robustsigmoid.png"),
       width = 10, height = 9, units="in", dpi=300)

################################################################################
# Figure B: pairwise raw distributions
################################################################################

# ASD raw data
ASD_pairwise_pyspi14_raw <- readRDS(paste0(ASD_rdata_path, "ABIDE_ASD_pyspi14_filtered.Rds"))
ASD_pairwise_pyspi14_raw %>%
  mutate(SPI = str_replace_all(SPI, "_", " ")) %>%
  ggplot(data=., mapping = aes(x=value, fill=SPI)) +
  geom_histogram() +
  ggtitle("ASD -- Raw") +
  facet_wrap(SPI ~ ., scales="free", ncol=4,
             labeller = labeller(SPI = label_wrap_gen(25))) +
  xlab("Feature Values") +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave(paste0(plot_path, "ASD_pyspi14_raw.png"),
       width = 9, height = 7.5, units="in", dpi=300)

# ASD z-scored data
ASD_pairwise_pyspi14_zscore <- readRDS(paste0(ASD_rdata_path, "ABIDE_ASD_pyspi14_filtered_zscored.Rds"))
ASD_pairwise_pyspi14_zscore %>%
  mutate(names = str_replace_all(names, "_", " ")) %>%
  ggplot(data=., mapping = aes(x=values, fill=names)) +
  geom_histogram() +
  ggtitle("ASD -- z-scored") +
  facet_wrap(names ~ ., scales="free", ncol=4,
             labeller = labeller(names = label_wrap_gen(25))) +
  xlab("Feature Values") +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave(paste0(plot_path, "ASD_pyspi14_zscored.png"),
       width = 9, height = 7.5, units="in", dpi=300)

# ASD robust sigmoid data
ASD_pairwise_pyspi14_robustsigmoid <- ASD_pairwise_pyspi14_raw %>%
  dplyr::rename("names" = "SPI",
                "values" = "value") %>%
  normalise_feature_frame(., method="RobustSigmoid") 

ASD_pairwise_pyspi14_robustsigmoid %>%
  mutate(names = str_replace_all(names, "_", " ")) %>%
  ggplot(data=., mapping = aes(x=values, fill=names)) +
  geom_histogram() +
  ggtitle("ASD -- Robust Sigmoid") +
  facet_wrap(names ~ ., scales="free", ncol=4,
             labeller = labeller(names = label_wrap_gen(25))) +
  xlab("Feature Values") +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave(paste0(plot_path, "ASD_pyspi14_robustsigmoid.png"),
       width = 9, height = 7.5, units="in", dpi=300)