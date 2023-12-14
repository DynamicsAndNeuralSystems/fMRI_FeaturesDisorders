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
catch25_feature_info <- read.csv("~/github/fMRI_FeaturesDisorders/data_visualisation/catch25_info.csv")
pyspi14_feature_info <- read.csv("~/github/fMRI_FeaturesDisorders/data_visualisation/SPI_info.csv")

mean_binned_in_R <-UCLA_CNP_catch25 %>%
  filter(names=="DN_Mean") %>%
  mutate(bin = cut_interval(values, n=100)) %>%
  group_by(names, bin) %>%
  count() %>%
  separate(bin, sep=",", into=c("lower_bound", "upper_bound")) %>%
  mutate(lower_bound = as.numeric(gsub("\\[|\\(", "", lower_bound)),
         upper_bound = as.numeric(gsub("\\]", "", upper_bound)))  %>%
  rowwise() %>%
  mutate(bin_center = mean(c(lower_bound, upper_bound)), .keep="unused") 

mean_binned_in_python = pyarrow_feather$read_feather("~/data/TS_feature_manuscript/all_normalisations_counts.feather") %>%
  filter(names=="DN_Mean",
         Normalisation=="Raw_Values") %>%
  separate(bin, sep=",", into=c("lower_bound", "upper_bound")) %>%
  mutate(lower_bound = as.numeric(gsub("\\[", "", lower_bound)),
         upper_bound = as.numeric(gsub("\\]", "", upper_bound))) %>%
  rowwise() %>%
  mutate(bin_center = mean(c(lower_bound, upper_bound)), .keep="unused") 

normalized_bin_data %>%
  ggplot(data=., mapping=aes(x=bin_center, y=count, fill=names)) +
  geom_bar(stat="identity")

# Load univariate catch25 data
UCLA_CNP_catch25_raw_data <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_catch25_filtered.feather")) %>%
  mutate(Normalization = "Raw Data")

ABIDE_ASD_catch25_raw_data <- pyarrow_feather$read_feather(glue("{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_catch25_filtered.feather")) %>%
  mutate(Normalization = "Raw Data")

UCLA_CNP_catch25_z <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_catch25_filtered_zscored.feather")) %>%
  mutate(Normalization = "z-scored")

ABIDE_ASD_catch25_z <- pyarrow_feather$read_feather(glue("{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_catch25_filtered_zscored.feather")) %>%
  mutate(Normalization = "z-score")

UCLA_CNP_catch25_RS <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_catch25_filtered_MixedSigmoid.feather")) %>%
  mutate(Normalization = "Mixed Sigmoid")

ABIDE_ASD_catch25_RS <- pyarrow_feather$read_feather(glue("{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_catch25_filtered_MixedSigmoid.feather")) %>%
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
  mutate(Normalization = "MixedSigmoid")

ABIDE_ASD_pyspi14_RS <- pyarrow_feather$read_feather(glue("{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_pyspi14_filtered_MixedSigmoid.feather")) %>%
  mutate(Normalization = "MixedSigmoid")


# Plot values for UCLA CNP
# Function to plot raw data, z-scored data, and mixed sigmoid-transformed data for each catch25 featuree
plot_values <- function(binned_feature_data, TS_feature_info, norm_type="Raw_Values", y_label="Raw\nValues") {
  
  binned_feature_data %>%
    group_by(feature_name) %>%
    mutate(feature_bin = cut(values_rounded, breaks = seq(0, max(values_rounded) + max(values_rounded)/100, max(values_rounded)/100))) %>%
    group_by(feature_name, feature_bin) %>%
    summarise(total_count = sum(count))
  
  p <- binned_feature_data %>%
    left_join(., TS_feature_info) %>%
    filter(Normalisation == norm_type,Figure_name=="Mean") %>%
    group_by(feature_name) %>%
    mutate(feature_bin = cut(values_rounded, breaks = seq(0, max(values_rounded) + max(values_rounded)/10, max(values_rounded)/10))) %>%
    group_by(feature_name, feature_bin) %>%
    summarise(total_count = sum(count)) %>%
    filter(total_count != 0, is.na(feature_bin)) %>%
    mutate(Figure_name = gsub("_", " ", Figure_name)) %>%
    ggplot(data = ., mapping = aes(x = feature_bin,
                                   y = total_count, 
                                   fill = Figure_name)) +
    geom_bar(stat="identity") +
    facet_wrap(Figure_name ~ ., scales="free", nrow=1,
               labeller = labeller(Figure_name = label_wrap_gen(10))) +
    ylab(y_label) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=10, angle=45, hjust=1),
          # axis.text.y = element_blank(),
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

# catch25
UCLA_CNP_catch25_plot_list <- list(plot_values(feature_data=UCLA_CNP_catch25_raw_data %>% dplyr::rename("feature_name" = "names"), TS_feature_info=catch25_feature_info, norm_type="none", "Raw\nValues"),
                                   plot_values(UCLA_CNP_catch25_z %>% dplyr::rename("feature_name" = "names"), TS_feature_info=catch25_feature_info, norm_type="z_score", "z-score"),
                                   plot_values(UCLA_CNP_catch25_RS %>% dplyr::rename("feature_name" = "names"), TS_feature_info=catch25_feature_info, norm_type="mixed_sigmoid", "Mixed\nSigmoid"))

wrap_plots(UCLA_CNP_catch25_plot_list, ncol=1, heights=c(0.36, 0.3, 0.3))
ggsave(glue("{plot_path}/UCLA_CNP_catch25_norms.svg"), bg="white",
       width = 28, height = 5, units = "in", dpi = 300)

# pyspi14
UCLA_CNP_pyspi14_plot_list <- list(plot_values(UCLA_CNP_pyspi14_raw_data %>% dplyr::rename("pyspi_name" = "names"), TS_feature_info=pyspi14_feature_info, norm_type="none", "Raw\nValues"),
                                   plot_values(UCLA_CNP_pyspi14_z %>% dplyr::rename("pyspi_name" = "names"), TS_feature_info=pyspi14_feature_info, norm_type="z_score", "z-score"),
                                   plot_values(UCLA_CNP_pyspi14_RS %>% dplyr::rename("pyspi_name" = "names"), TS_feature_info=pyspi14_feature_info, norm_type="mixed_sigmoid", "Mixed\nSigmoid"))

wrap_plots(UCLA_CNP_pyspi14_plot_list, ncol=1, heights=c(0.36, 0.3, 0.3))
ggsave(glue("{plot_path}/UCLA_CNP_pyspi14_norms.svg"), bg="white",
       width = 28, height = 5, units = "in", dpi = 300)

# Plot values for ABIDE ASD

# catch25
ABIDE_ASD_catch25_plot_list <- list(plot_values(ABIDE_ASD_catch25_raw_data %>% dplyr::rename("feature_name" = "names"), TS_feature_info=catch25_feature_info, norm_type="none", "Raw\nValues"),
                                    plot_values(ABIDE_ASD_catch25_z %>% dplyr::rename("feature_name" = "names"), TS_feature_info=catch25_feature_info, norm_type="z_score", "z-score"),
                                    plot_values(ABIDE_ASD_catch25_RS %>% dplyr::rename("feature_name" = "names"), TS_feature_info=catch25_feature_info, norm_type="mixed_sigmoid", "Mixed\nSigmoid"))

wrap_plots(ABIDE_ASD_catch25_plot_list, ncol=1, heights=c(0.36, 0.3, 0.3))
ggsave(glue("{plot_path}/ABIDE_ASD_catch25_norms.svg"), bg="white",
       width = 28, height = 5, units = "in", dpi = 300)

# pyspi14
ABIDE_ASD_pyspi14_plot_list <- list(plot_values(ABIDE_ASD_pyspi14_raw_data %>% dplyr::rename("pyspi_name" = "names"), TS_feature_info=pyspi14_feature_info, norm_type="none", "Raw\nValues"),
                                    plot_values(ABIDE_ASD_pyspi14_z %>% dplyr::rename("pyspi_name" = "names"), TS_feature_info=pyspi14_feature_info, norm_type="z_score", "z-score"),
                                    plot_values(ABIDE_ASD_pyspi14_RS %>% dplyr::rename("pyspi_name" = "names"), TS_feature_info=pyspi14_feature_info, norm_type="mixed_sigmoid", "Mixed\nSigmoid"))

wrap_plots(ABIDE_ASD_pyspi14_plot_list, ncol=1, heights=c(0.36, 0.3, 0.3))
ggsave(glue("{plot_path}/ABIDE_ASD_pyspi14_norms.svg"), bg="white",
       width = 28, height = 5, units = "in", dpi = 300)