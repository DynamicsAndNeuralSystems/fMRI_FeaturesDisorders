python_to_use <- "~/.conda/envs/pyspi/bin/python3"

reticulate::use_python(python_to_use)
library(feather)
library(reticulate)
library(tidyverse)
library(theft)
library(glue)
library(patchwork)
library(cowplot)
theme_set(theme_cowplot())

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

# Import normalisation methods from python
source_python("~/github/fMRI_FeaturesDisorders/helper_functions/classification/mixed_sigmoid_normalisation.py")

# DIY rlist::list.append
list.append <- function (.data, ...) 
{
  if (is.list(.data)) {
    c(.data, list(...))
  }
  else {
    c(.data, ..., recursive = FALSE)
  }
}

# Define plot path
plot_path <- "~/github/fMRI_FeaturesDisorders/plots/Manuscript_Draft/normalisation_analysis/"
TAF::mkdir(plot_path)

# Define data paths
UCLA_CNP_data_path <- "~/data/UCLA_CNP/processed_data/"
ABIDE_ASD_data_path <- "~/data/ABIDE_ASD/processed_data/"

# Load feature info data
TS_feature_info <- read.csv(glue("{github_dir}/data_visualisation/catch24_info.csv"))

# Load univariate catch24 data
UCLA_CNP_catch24_raw_data <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_catch24_filtered.feather")) %>%
  mutate(Normalization = "Raw Data")

ABIDE_ASD_catch24_raw_data <- pyarrow_feather$read_feather(glue("{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_catch24_filtered.feather")) %>%
  mutate(Normalization = "Raw Data")

UCLA_CNP_catch24_z <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_catch24_filtered_zscored.feather")) %>%
  mutate(Normalization = "z-scored")

ABIDE_ASD_catch24_z <- pyarrow_feather$read_feather(glue("{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_catch24_filtered_zscored.feather")) %>%
  mutate(Normalization = "z-score")

UCLA_CNP_catch24_RS <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_catch24_filtered_MixedSigmoid.feather")) %>%
  mutate(Normalization = "Mixed Sigmoid")

ABIDE_ASD_catch24_RS <- pyarrow_feather$read_feather(glue("{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_catch24_filtered_MixedSigmoid.feather")) %>%
  mutate(Normalization = "Mixed Sigmoid")

# Load pairwise pyspi14 data
UCLA_CNP_pyspi14_raw_data <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered.feather")) %>%
  mutate(Normalization = "Raw Data") %>%
  dplyr::rename("names"="SPI", "values"="value") %>%
  mutate(Brain_Region = paste0(brain_region_from, "_", brain_region_to), .keep = "unused")

ABIDE_ASD_pyspi14_raw_data <- pyarrow_feather$read_feather(glue("{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_pyspi14_filtered.feather")) %>%
  mutate(Normalization = "Raw Data") %>%
  dplyr::rename("names"="SPI", "values"="value") %>%
  mutate(Brain_Region = paste0(brain_region_from, "_", brain_region_to), .keep = "unused")

UCLA_CNP_pyspi14_z <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered_zscored.feather")) %>%
  mutate(Normalization = "z-scored")

ABIDE_ASD_pyspi14_z <- pyarrow_feather$read_feather(glue("{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_pyspi14_filtered_zscored.feather")) %>%
  mutate(Normalization = "z-score")

UCLA_CNP_pyspi14_RS <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered_MixedSigmoid.feather")) %>%
  mutate(Normalization = "Mixed Sigmoid")

ABIDE_ASD_pyspi14_RS <- pyarrow_feather$read_feather(glue("{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_pyspi14_filtered_MixedSigmoid.feather")) %>%
  mutate(Normalization = "Mixed Sigmoid")


# Plot values for UCLA CNP
# Function to plot raw data, z-scored data, and mixed sigmoid-transformed data for each catch24 featuree
plot_values <- function(feature_data, norm_type="none", y_label="Raw\nValues") {
  p <- feature_data %>%
    dplyr::rename("TS_Feature" = "names") %>%
    left_join(., TS_feature_info) %>%
    ggplot(data = ., mapping = aes(x = values,
                                   y = after_stat(count)/sum(after_stat(count)), 
                                   fill = Figure_Name)) +
    geom_histogram() +
    facet_wrap(Figure_Name ~ ., scales="free", nrow=1,
               labeller = labeller(Figure_Name = label_wrap_gen(10))) +
    ylab(y_label) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=10, angle=45, hjust=1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_text(angle=0, vjust=0.5, size=16))
  
  if (norm_type != "none") {
    p <- p + 
      theme(strip.background = element_blank(),
            strip.text = element_blank())
  }
  
  if (norm_type == "mixed_sigmoid") {
    p <- p + 
      scale_x_continuous(breaks = c(0, 0.5, 1))
  }
  return(p)
}

# catch24
UCLA_CNP_catch24_plot_list <- list(plot_values(UCLA_CNP_catch24_raw_data, norm_type="none", "Raw\nValues"),
                                   plot_values(UCLA_CNP_catch24_z, norm_type="z_score", "z-score"),
                                   plot_values(UCLA_CNP_catch24_RS, norm_type="mixed_sigmoid", "Mixed\nSigmoid"))

wrap_plots(UCLA_CNP_catch24_plot_list, ncol=1, heights=c(0.36, 0.3, 0.3))
ggsave(glue("{plot_path}/UCLA_CNP_catch24_norms.png"), bg="white",
       width = 28, height = 5, units = "in", dpi = 300)

# pyspi14
UCLA_CNP_pyspi14_plot_list <- list(plot_values(UCLA_CNP_pyspi14_raw_data, norm_type="none", "Raw\nValues"),
                                   plot_values(UCLA_CNP_pyspi14_z, norm_type="z_score", "z-score"),
                                   plot_values(UCLA_CNP_pyspi14_RS, norm_type="mixed_sigmoid", "Mixed\nSigmoid"))

wrap_plots(UCLA_CNP_pyspi14_plot_list, ncol=1, heights=c(0.36, 0.3, 0.3))
ggsave(glue("{plot_path}/UCLA_CNP_pyspi14_norms.png"), bg="white",
       width = 28, height = 5, units = "in", dpi = 300)

# Plot values for ABIDE ASD

# catch24
ABIDE_ASD_catch24_plot_list <- list(plot_values(ABIDE_ASD_catch24_raw_data, norm_type="none", "Raw\nValues"),
                                    plot_values(ABIDE_ASD_catch24_z, norm_type="z_score", "z-score"),
                                    plot_values(ABIDE_ASD_catch24_RS, norm_type="mixed_sigmoid", "Mixed\nSigmoid"))

wrap_plots(ABIDE_ASD_catch24_plot_list, ncol=1, heights=c(0.36, 0.3, 0.3))
ggsave(glue("{plot_path}/ABIDE_ASD_catch24_norms.png"), bg="white",
       width = 28, height = 5, units = "in", dpi = 300)

# pyspi14
ABIDE_ASD_pyspi14_plot_list <- list(plot_values(ABIDE_ASD_pyspi14_raw_data, norm_type="none", "Raw\nValues"),
                                    plot_values(ABIDE_ASD_pyspi14_z, norm_type="z_score", "z-score"),
                                    plot_values(ABIDE_ASD_pyspi14_RS, norm_type="mixed_sigmoid", "Mixed\nSigmoid"))

wrap_plots(ABIDE_ASD_pyspi14_plot_list, ncol=1, heights=c(0.36, 0.3, 0.3))
ggsave(glue("{plot_path}/ABIDE_ASD_pyspi14_norms.png"), bg="white",
       width = 28, height = 5, units = "in", dpi = 300)