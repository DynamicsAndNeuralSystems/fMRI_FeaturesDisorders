#------------------------------------
# This script sets out to produce a
# function for reading in matlab time
# series files into R
#------------------------------------

#--------------------------------------
# Author: Trent Henderson, 9 March 2021
# Updated: Annie Bryant, 15 March 2022
#--------------------------------------

require(plyr)
library(tidyverse)
library(rmatio)
library(R.matlab)
library(theft)
library(rlist)
library(argparse)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Source scripts
setwd(getSrcDirectory()[1])
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("feature_calculation/catch22_all_regions.R")
source("feature_analysis/region_by_region_analysis.R")
source("visualization/plot_feature_acc_across_ROIs.R")
source("visualization/violin_plots.R")
source("prep_data/load_mat_data.R")
source("prep_data/get_dx_breakdown.R")
source("prep_data/compile_movement_data.R")
set.seed(127)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Parse command-line arguments
parser <- ArgumentParser(description='Read in matlab time-series data and convert to an R data object.')
parser$add_argument("--mat_file", help=".mat file containing the time-series data and other metadata.")
parser$add_argument("--label_metadata", help="CSV file containing sample metadata info.",
                    default="D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/participants.csv")
parser$add_argument("--data_path", help="File path to store resulting Rdata objects.",
                    default="D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/")
parser$add_argument("--plot_path", help="File path to store plot images.",
                    default="D:/Virtual_Machines/Shared_Folder/PhD_work/plots/")
parser$add_argument("--overwrite", help="Should the Rdata object be overwritten if it already exists? Default is F.",
                    action="store_true", default=FALSE)

# Parse arguments
args <- parser$parse_args()
mat_file <- args$mat_file
label_metadata <- args$label_metadata
data_path <- args$data_path
plot_path <- args$plot_path
overwrite <- args$overwrite
rdata_path <- paste0(data_path, "Rdata/")

# DEBUG ONLY
# data_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/"
# rdata_path <- paste0(data_path, "Rdata/")
# mat_file <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/new/UCLA_time_series_four_groups.mat"
# label_metadata <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/participants.csv"
# subject_csv <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/participants.csv"
# plot_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/plots/"
# overwrite <- T
# UCLA_AROMA_2P_catch22_ROIwise_t_test <- readRDS(paste0(rdata_path, "UCLA_AROMA_2P_catch22_ROIwise_t_test.Rds"))
# UCLA_AROMA_2P_catch22 <- readRDS(paste0(rdata_path, "UCLA_AROMA_2P_catch22.Rds"))
# 
# 


# Plot top 8 features for an example ROI
violin_plot_for_ROI(feature_matrix = feature_matrix,
                        class_res = class_res,
                        this_ROI = this_ROI,
                        num_feature = num_feature)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Plot the distribution of t-scores across features by noise processing method
for (noise_proc in c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")) { 
  
  # Clean up names
  noise_label <- gsub("\\+", "_", noise_proc)# Load catch22 feature matrix
  
  class_res <- readRDS(paste0(rdata_path, "UCLA_", 
                              noise_label, "_catch22_ROIwise_",  
                              test_label, ".Rds"))
  
  plot_ROI_acc_by_feature(region_wise_univ_class_res = class_res,
                          xlab = "T statistic",
                          noise_proc = noise_label,
                          plot_path = plot_path)
}


#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Feature-by-feature analysis
this_feature <- "DN_HistogramMode_5"
feature_matrix <- UCLA_AROMA_2P_catch22
class_res <- UCLA_AROMA_2P_catch22_ROIwise_t_test
num_ROI <- 8

PCA_dimplot_for_ROI(feature_matrix = feature_matrix,
                    class_res = class_res,
                    this_ROI = this_ROI)
ggsave(paste0(plot_path, "PC1_PC2_left_bankssts.png"),
       width = 6, height = 5, units="in", dpi=300)

violin_plot_for_feature(feature_matrix = feature_matrix,
                        class_res = class_res,
                        this_feature = this_feature,
                        num_ROI = num_ROI)

