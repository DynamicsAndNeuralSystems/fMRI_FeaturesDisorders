#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------

# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--github_dir", default="/headnode1/abry4213/github/")
parser$add_argument("--data_path", default="/headnode1/abry4213/data/UCLA_Schizophrenia/")
parser$add_argument("--pkl_file", default="calc.pkl")
parser$add_argument("--python_to_use", default="/headnode1/abry4213/.conda/envs/pyspi/bin/python3")
parser$add_argument("--univariate_feature_set", default="catch22")
parser$add_argument("--pairwise_feature_set", default="pyspi14")
parser$add_argument("--sample_metadata_file", default="UCLA_Schizophrenia_sample_metadata.Rds")
parser$add_argument("--brain_region_lookup", default="", nargs='?')
parser$add_argument("--noise_procs", default=c(""))
parser$add_argument("--main_noise_proc", default="AROMA+2P+GMR")
parser$add_argument("--dataset_ID", default="UCLA_Schizophrenia")

# Parse input arguments
args <- parser$parse_args()
github_dir <- args$github_dir
data_path <- args$data_path
pkl_file <- args$pkl_file
python_to_use <- args$python_to_use
univariate_feature_set <- args$univariate_feature_set
pairwise_feature_set <- args$pairwise_feature_set
sample_metadata_file <- args$sample_metadata_file
brain_region_lookup <- args$brain_region_lookup
noise_procs <- args$noise_procs
main_noise_proc <- args$main_noise_proc
dataset_ID <- args$dataset_ID

# python_to_use <- "/headnode1/abry4213/.conda/envs/pyspi/bin/python3"
# univariate_feature_set <- "catch22"
# pairwise_feature_set <- "pyspi14_mod"
# github_dir <- "/headnode1/abry4213/github/"
# pkl_file <- "calc.pkl"

# UCLA schizophrenia
# data_path <- "/headnode1/abry4213/data/UCLA_Schizophrenia/"
# sample_metadata_file <- "UCLA_Schizophrenia_sample_metadata.Rds"
# dataset_ID <- "UCLA_Schizophrenia"
# noise_procs <- c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")
# brain_region_lookup <- "Brain_Region_info.csv"

# ABIDE ASD
# data_path <- "/headnode1/abry4213/data/ABIDE_ASD/"
# sample_metadata_file <- "ABIDE_ASD_sample_metadata.Rds"
# dataset_ID <- "ABIDE_ASD"
# noise_procs <- c("FC1000")
# brain_region_lookup <- "Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv"

# HCP100
# data_path <- "/headnode1/abry4213/data/HCP100/"
# sample_metadata_file <- "HCP100_sample_metadata.Rds"
# dataset_ID <- "HCP100"
# noise_procs <- c("AROMA+2P+GMR")
# brain_region_lookup <- "Brain_Region_info.csv"

sample_metadata <- readRDS(paste0(data_path, sample_metadata_file))
pkl_data_path <- paste0(data_path, "raw_data/numpy_files/")

rdata_path <- paste0(data_path, "processed_data/Rdata/")
plot_dir <- paste0(data_path, "plots/")

icesTAF::mkdir(rdata_path)
icesTAF::mkdir(plot_dir)

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

# Load specified python version
reticulate::use_python(python_to_use)
library(reticulate)
reticulate::source_python(paste0(github_dir, "fMRI_FeaturesDisorders/helper_functions/data_prep_and_QC/pickle_reader.py"))
library(tidyverse)

# Source QC functions
source(paste0(github_dir, "fMRI_FeaturesDisorders/helper_functions/data_prep_and_QC/QC_functions_pairwise.R"))

# Unlist noise-processing methods
tryCatch({
  noise_procs <- stringr::str_split(noise_procs, ";")[[1]]
  noise_procs <- unlist(noise_procs)
}, error = function(e) {
  message(e)
})
# Remove empty noise-processing methods
noise_procs <- noise_procs[noise_procs!=""]

# Print out arguments
cat("pkl_data_path:", pkl_data_path, "\n")
cat("rdata_path:", rdata_path, "\n")
cat("noise_procs:", paste(noise_procs, collapse=", "), "\n")

#-------------------------------------------------------------------------------
# Function to read in pyspi pickle files per sample and convert
# The SPI result data into an RDS file per sample
#-------------------------------------------------------------------------------
read_pyspi_pkl_into_RDS <- function(pkl_data_path,
                                    rdata_path,
                                    sample_metadata,
                                    noise_procs = c("AROMA+2P",
                                                    "AROMA+2P+GMR",
                                                    "AROMA+2P+DiCER")) {
  # Iterate over each noise-processing method
  for (noise_proc in noise_procs) {
    noise_label = gsub("\\+", "_", noise_proc)
    cat("\nNow processing", noise_label, "\n")
    
    # Define data path for this noise processing method
    np_data_path <- paste0(pkl_data_path, noise_label, "/")
    
    # Define output rdata path
    np_rdata_path <- paste0(rdata_path, noise_label, "_", pairwise_feature_set, "/")
    icesTAF::mkdir(np_rdata_path)
    
    # Iterate over each sample
    for (sample in unique(list.dirs(np_data_path, recursive = F, full.names = F))) {
      
      # If sample doesn't have a corresponding pyspi RDS file for this 
      # noise-processing method, create one
      if (!file.exists(paste0(np_rdata_path, sample, "_pyspi.Rds"))) {
        cat("\nNow prepping data for", sample, noise_label, "\n")
        tryCatch({
          sample_pkl_data <- extract_df_from_pkl(paste0(np_data_path, sample, "/", pkl_file)) %>%
            mutate(Sample_ID = sample,
                   Noise_Proc = noise_proc,
                   brain_region_1 = as.numeric(gsub("proc-", "", brain_region_1)),
                   brain_region_2 = as.numeric(gsub("proc-", "", brain_region_2)))
          
          if ("Diagnosis" %in% colnames(sample_metadata)) {
            sample_pkl_data <- sample_pkl_data %>%
              mutate(Diagnosis = subset(sample_metadata, Sample_ID == sample) %>% pull(Diagnosis))
          }
          # Save results to an RDS file for this sample
          saveRDS(sample_pkl_data, file=paste0(np_rdata_path, sample, "_pyspi.Rds"))},
          error = function(e){
            cat("\nError for sample", sample, "\n")
            print(e)
          })
      }
    }
  }
}
#-------------------------------------------------------------------------------
# Function to merge all of the individual sample pyspi .Rds files into one
#-------------------------------------------------------------------------------
merge_pyspi_res_for_study <- function(rdata_path,
                                      dataset_ID = "UCLA_Schizophrenia",
                                      brain_region_lookup,
                                      noise_procs = c("AROMA+2P",
                                                      "AROMA+2P+GMR",
                                                      "AROMA+2P+DiCER")) {
  
  # Read in ROI index data
  ROI_index <- read.csv(brain_region_lookup)
  
  if (!file.exists(paste0(rdata_path, dataset_ID, 
                          "_", pairwise_feature_set, ".Rds"))) {
    # Create list to store results across noise-processing methods
    noise_proc_res <- list()
    # Iterate over each noise-processing method
    for (noise_proc in noise_procs) {
      
      noise_label = gsub("\\+", "_", noise_proc)
      # Define data path for this noise processing method
      np_rdata_path <- paste0(rdata_path, noise_label, "_", pairwise_feature_set, "/")
      
      # Iterate over each sample
      pyspi_files = unique(list.files(np_rdata_path, recursive = F, full.names = F))
      for (file in pyspi_files) {
        cat("Now processing", file, "\n")
        sample = gsub(sprintf("_%s|.Rds", pairwise_feature_set), "", file)
        
        # If sample doesn't have a corresponding pyspi RDS file for this 
        # noise-processing method, create one
        tryCatch({
          sample_pyspi_res <- readRDS(paste0(np_rdata_path, file)) 
          
          # Append results to list
          noise_proc_res <- list.append(noise_proc_res, sample_pyspi_res)
        },
        error = function(e){
          cat("\nError for sample", sample, "\n")
          print(e)
        })
      }
    }
    
    # Merge subjects' data together
    all_pyspi_data <- do.call(plyr::rbind.fill, noise_proc_res)  %>%
      mutate(comparison = row_number(),
             group = stringr::str_to_sentence(Diagnosis)) %>%
      pivot_longer(cols = c(brain_region_1,
                            brain_region_2),
                   names_to = "Region_Number",
                   values_to = "Index") %>%
      # Convert Index to number and add 1 since python is base 0
      # While R is base 1
      mutate(Index = 1 + as.numeric(gsub("proc-", "", Index))) %>%
      left_join(ROI_index) %>%
      dplyr::select(-Index) %>%
      pivot_wider(id_cols = c("Sample_ID", "Diagnosis", "SPI", "Noise_Proc", "value", "comparison"),
                  names_from = "Region_Number",
                  values_from = "Brain_Region") %>%
      dplyr::select(-comparison) 
    
    saveRDS(all_pyspi_data, paste0(rdata_path, dataset_ID, 
                                   "_", pairwise_feature_set, ".Rds"))
    
  }
}

#-------------------------------------------------------------------------------
# Call functions
#-------------------------------------------------------------------------------
read_pyspi_pkl_into_RDS(pkl_data_path = pkl_data_path,
                        rdata_path = rdata_path,
                        sample_metadata = sample_metadata,
                        noise_procs = noise_procs)

merge_pyspi_res_for_study(rdata_path = rdata_path,
                          dataset_ID = dataset_ID,
                          brain_region_lookup = paste0(data_path, brain_region_lookup),
                          noise_procs = noise_procs)

#-------------------------------------------------------------------------------
# Perform QC for pairwise data
#-------------------------------------------------------------------------------
run_QC_for_dataset(data_path = data_path, 
                   proc_rdata_path = rdata_path,
                   sample_metadata_file = sample_metadata_file,
                   dataset_ID = dataset_ID,
                   pairwise_feature_set = pairwise_feature_set,
                   noise_proc = main_noise_proc)