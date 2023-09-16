#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------

univariate_feature_set <- "catch24"
github_dir <- "~/Library/CloudStorage/OneDrive-TheUniversityofSydney(Students)/github/"

# UCLA CNP
data_path <- "~/data/UCLA_CNP/"
dataset_ID <- "UCLA_CNP"
sample_metadata_file <- "UCLA_CNP_sample_metadata.feather"
noise_procs <- c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")
brain_region_lookup <- "UCLA_CNP_Brain_Region_Lookup.feather"

rdata_path <- paste0(data_path, "processed_data/")
plot_dir <- paste0(data_path, "plots/")

TAF::mkdir(plot_dir)
TAF::mkdir(rdata_path)

# Set the seed
set.seed(127)

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"
reticulate::use_python(python_to_use)

# Load tidyverse
library(tidyverse)
library(theft)
library(feather)
library(reticulate)
library(cowplot)
theme_set(theme_cowplot())

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

#-------------------------------------------------------------------------------
# Source helper scripts
#-------------------------------------------------------------------------------
# Set working directory to file location
helper_script_dir = paste0(github_dir, "fMRI_FeaturesDisorders/helper_functions/")
source(paste0(helper_script_dir, "data_prep_and_QC/TS_feature_extraction.R"))
source(paste0(helper_script_dir, "data_prep_and_QC/QC_functions_univariate.R"))

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



#-------------------------------------------------------------------------------
# Read in TS data for dataset, partitioned by noise-processing method
#-------------------------------------------------------------------------------

# Load brain region lookup table

brain_region_lookup_table <- pyarrow_feather$read_feather(paste0(data_path, "study_metadata/", brain_region_lookup))

read_in_sample_TS_data <- function(sample_ID, dataset_ID, noise_proc,
                                   brain_region_lookup_table) {
  noise_label <- gsub("\\+", "_", noise_proc)
  tryCatch({TS_data <- read.csv(paste0(data_path,
                                       "raw_data/time_series_files/",
                                       noise_label, "/",
                                       sample_ID, "_TS.csv"),
                                header=T) %>%
    mutate(timepoint = 1:nrow(.)) %>%
    pivot_longer(cols = c(-timepoint),
                 names_to = "Index",
                 values_to = "values") %>%
    mutate(Sample_ID = sample_ID,
           Noise_Proc = noise_proc,
           Index = as.numeric(gsub("X|V", "", Index))) %>%
    left_join(., brain_region_lookup_table) %>%
    dplyr::select(Sample_ID, Noise_Proc, Brain_Region, timepoint, values)
  return(TS_data)},
  error = function(e) {
    cat("\nError for subject", sample_ID, "\n")
    message(e)
  })
}

noise_proc_res_list <- list()

# Iterate over each of the three noise proc methods
for (noise_proc in noise_procs) {
  noise_label <- gsub("\\+", "_", noise_proc)
    # filename-friendly noise proc label
    noise_label <- gsub("\\+", "_", noise_proc)
    # Find sample IDs
    sample_IDs <- list.files(paste0(data_path, "raw_data/time_series_files/", noise_label)) %>%
      gsub("_TS.csv", "", .)
    # Read in time-series data across all participants
    np_TS_data <- sample_IDs %>%
      purrr::map_df(~ read_in_sample_TS_data(sample_ID = .x,
                                             dataset_ID = dataset_ID,
                                             noise_proc = noise_proc,
                                             brain_region_lookup_table = brain_region_lookup_table))
    
    # Append results to list
    noise_proc_res_list <- list.append(noise_proc_res_list, np_TS_data)
    
}

# Concatenate results
noise_proc_res <- do.call(plyr::rbind.fill, noise_proc_res_list)

#-------------------------------------------------------------------------------
# Check for any subjects with flat time series across all regions
#-------------------------------------------------------------------------------
all_NA_samples <- noise_proc_res %>% 
  group_by(Sample_ID, Noise_Proc) %>%
  summarise(SD_values = sd(values, na.rm=T)) %>%
  filter(SD_values==0) %>% 
  distinct(Sample_ID) %>%
  pull(Sample_ID)

noise_proc_res %>%
  filter(Sample_ID %in% all_NA_samples) %>%
  ggplot(data=., mapping=aes(x=timepoint, y=values, color=Brain_Region)) +
  ggtitle(sprintf("Raw time-series for %s\nNA samples with %s",
                  gsub("_", " ", dataset_ID), univariate_feature_set)) +
  geom_line(alpha=0.6) +
  facet_grid(Sample_ID ~ Noise_Proc, switch="y") +
  theme(legend.position="none",
        strip.text.y.left = element_text(angle=0),
        plot.title = element_text(hjust=0.5))