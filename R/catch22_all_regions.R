#------------------------------------
# This script runs catch22 for all brain region across 
# schizophrenia vs. control subjects
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 15 March 2022
#--------------------------------------

require(plyr)
library(tidyverse)
library(theft)
library(argparse)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
parser <- ArgumentParser(description='Runs catch22 for one brain region across all schizophrenia vs. control subjects.')
parser$add_argument("--rdata_path", help="File path containing Rdata objects with time-series data.",
                    default="D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/Rdata/")
parser$add_argument("--noise_proc", help="Noise processing method from which data should be extracted.")

# Parse arguments
args <- parser$parse_args()
rdata_path <- args$rdata_path
noise_proc <- args$noise_proc
noise_label <- gsub("\\+", "_", noise_proc)

# DEBUG ONLY
rdata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/Rdata/"
noise_proc <- "AROMA+2P"
noise_label <- gsub("\\+", "_", noise_proc)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

# Load time-series data
TS_data <- readRDS(paste0(rdata_path, sprintf("UCLA_%s.rds", noise_label)))

# Create a new ID that contains both subject ID and brain region
TS_data$unique_ID <- paste(TS_data$Subject_ID, TS_data$Brain_Region, sep="_")

# Run catch22 using theft
feature_matrix <- calculate_features(data = TS_data, 
                                     id_var = "unique_ID", 
                                     time_var = "timepoint", 
                                     values_var = "value", 
                                     group_var = "diagnosis", 
                                     feature_set = "catch22",
                                     catch24 = F)