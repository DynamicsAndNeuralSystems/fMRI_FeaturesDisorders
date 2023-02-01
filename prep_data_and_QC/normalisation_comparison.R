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
source_python("~/github/fMRI_FeaturesDisorders/helper_functions/classification/robust_sigmoid_normalisation.py")

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
plot_path <- "~/github/fMRI_FeaturesDisorders/plots/Manuscript_Draft/FigureS2"
TAF::mkdir(plot_path)

# Define data paths
UCLA_CNP_data_path <- "~/data/UCLA_CNP/processed_data/"
ABIDE_ASD_data_path <- "~/data/ABIDE_ASD/processed_data/"

# Load metadata
UCLA_CNP_metadata <- pyarrow_feather$read_feather("~/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather")
ABIDE_ASD_metadata <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_sample_metadata.feather")

# Load catch22 data
UCLA_CNP_catch22_data <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/UCLA_CNP_AROMA_2P_GMR_catch22_filtered.feather")) %>%
  mutate(Normalization = "Raw Data")
ABIDE_ASD_catch22_data <- pyarrow_feather$read_feather(glue("{ABIDE_ASD_data_path}/ABIDE_ASD_FC1000_catch22_filtered.feather")) %>%
  mutate(Normalization = "Raw Data")

# Apply normalization
apply_transform_to_dataset <- function(dataset, metadata, transform_name) {
  # Initialise list for transformed data results
  transformed_data_list <- list()
  
  # Combine with metadata
  dataset <- left_join(dataset, metadata)
  
  # Iterate over each diagnostic group
  for (dx_group in unique(dataset$Diagnosis)) {
    dx_group_data <- subset(dataset, Diagnosis==dx_group)
    
    # Iterate over each brain region
    for (brain_region in unique(dx_group_data$Brain_Region)) {
      region_data <- dx_group_data %>%
        filter(Brain_Region == brain_region) %>%
        dplyr::select(Sample_ID, names, values) %>%
        tidyr::pivot_wider(id_cols = Sample_ID, names_from = names, values_from = values)
      
      subjects <- region_data$Sample_ID
      
      region_matrix <- region_data %>%
        dplyr::select(-Sample_ID) %>%
        as.matrix()
      
      TS_features <- colnames(region_matrix)
      
      # Instantiate the transformer
      if (transform_name == "z-score") {
        transformer <- StandardScaler()$fit(region_matrix)
      } else {
        transformer <- RobustSigmoidScaler(unit_variance=TRUE)$fit(region_matrix)
      }
      
      # Apply transformer
      region_data_trans <- as.data.frame(transformer$transform(region_matrix))
      colnames(region_data_trans) <- TS_features
      region_data_trans$Sample_ID <- subjects
      
      # Reshape from wide to long
      region_data_trans_long <- region_data_trans %>%
        pivot_longer(cols=c(-Sample_ID), names_to="names", values_to="values") %>%
        mutate(Brain_Region = brain_region,
               Diagnosis = dx_group)
      
      # Append region data to list
      transformed_data_list <- list.append(transformed_data_list, 
                                           region_data_trans_long)
    }
  }
  
  
  # Concatenate results along list
  transformed_data <- do.call(plyr::rbind.fill, transformed_data_list)
  
  return(transformed_data)
}

# z-score
UCLA_CNP_catch22_data_zscore <- apply_transform_to_dataset(UCLA_CNP_catch22_data,
                                                           UCLA_CNP_metadata,
                                                           "z-score") %>%
  mutate(Normalization = "z-score")

ABIDE_ASD_catch22_data_zscore <- apply_transform_to_dataset(ABIDE_ASD_catch22_data,
                                                            ABIDE_ASD_metadata,
                                                            "z-score") %>%
  mutate(Normalization = "z-score")

# robust sigmoid
UCLA_CNP_catch22_data_robust_sigmoid <- apply_transform_to_dataset(UCLA_CNP_catch22_data,
                                                                   UCLA_CNP_metadata,
                                                                   "RobustSigmoid") %>%
  mutate(Normalization = "Robust Sigmoid")
ABIDE_ASD_catch22_data_robust_sigmoid <- apply_transform_to_dataset(ABIDE_ASD_catch22_data,
                                                                    ABIDE_ASD_metadata,
                                                                    "RobustSigmoid") %>%
  mutate(Normalization = "Robust Sigmoid")

# Function to plot raw data, z-scored data, and robust sigmoid-transformed data for each catch22 featuree
plot_values <- function(feature_data, dx_group, metadata, norm_type="none", y_label="Raw\nValues") {
  p <- feature_data %>%
    left_join(., metadata) %>%
    filter(Diagnosis == dx_group) %>%
    mutate(names = gsub("_", " ", names)) %>%
    ggplot(data = ., mapping = aes(x = values,
                                   y = after_stat(count)/sum(after_stat(count)), 
                                   fill = names)) +
    geom_histogram() +
    facet_wrap(names ~ ., scales="free", nrow=1,
               labeller = labeller(names = label_wrap_gen(10))) +
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
  
  if (norm_type == "robust_sigmoid") {
    p <- p + 
      scale_x_continuous(breaks = c(0, 0.5, 1))
  }
  return(p)
}

# Plot values for UCLA CNP
for (dx in c("Schizophrenia", "ADHD", "Bipolar", "Control")) {
  UCLA_CNP_dx_plot_list <- list(plot_values(UCLA_CNP_catch22_data, dx_group = dx, metadata = UCLA_CNP_metadata, norm_type="none", "Raw\nValues"),
                                plot_values(UCLA_CNP_catch22_data_zscore, dx_group = dx, metadata = UCLA_CNP_metadata, norm_type="z_score", "z-score"),
                                plot_values(UCLA_CNP_catch22_data_robust_sigmoid, dx_group = dx, metadata = UCLA_CNP_metadata, norm_type="robust_sigmoid", "Robust\nSigmoid"))
  
  wrap_plots(UCLA_CNP_dx_plot_list, ncol=1, heights=c(0.36, 0.3, 0.3))
  ggsave(glue("{plot_path}/UCLA_CNP_catch22_norms_{dx}.png"), bg="white",
         width = 28, height = 5, units = "in", dpi = 300)
}


# Plot values for ABIDE ASD
for (dx in c("ASD", "Control")) {
  ABIDE_ASD_dx_plot_list <- list(plot_values(ABIDE_ASD_catch22_data, dx_group = dx, metadata = ABIDE_ASD_metadata, norm_type="none", "Raw\nValues"),
                              plot_values(ABIDE_ASD_catch22_data_zscore, dx_group = dx, metadata = ABIDE_ASD_metadata, norm_type="z_score", "z-score"),
                              plot_values(ABIDE_ASD_catch22_data_robust_sigmoid, dx_group = dx, metadata = ABIDE_ASD_metadata, norm_type="robust_sigmoid", "Robust\nSigmoid"))
  
  wrap_plots(ABIDE_ASD_dx_plot_list, ncol=1, heights=c(0.36, 0.3, 0.3))
  ggsave(glue("{plot_path}/ABIDE_ASD_catch22_norms.png"), bg="white",
         width = 28, height = 5, units = "in", dpi = 300)
}