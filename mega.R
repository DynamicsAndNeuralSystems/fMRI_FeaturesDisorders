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
source("prep_data/load_mat_data.R")
source("prep_data/get_dx_breakdown.R")
set.seed(127)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Parse command-line arguments
parser <- ArgumentParser(description='Read in matlab time-series data and convert to an R data object.')
parser$add_argument("--mat_file", help=".mat file containing the time-series data and other metadata.")
parser$add_argument("--label_metadata", help="CSV file containing sample metadata info.",
                    default="D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/participants.csv")
parser$add_argument("--rdata_path", help="File path to store resulting Rdata objects.",
                    default="D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/Rdata/")
parser$add_argument("--plot_path", help="File path to store plot images.",
                    default="D:/Virtual_Machines/Shared_Folder/PhD_work/plots/")
parser$add_argument("--overwrite", help="Should the Rdata object be overwritten if it already exists? Default is F.",
                    action="store_true", default=FALSE)

# Parse arguments
args <- parser$parse_args()
mat_file <- args$mat_file
label_metadata <- args$label_metadata
rdata_path <- args$rdata_path
plot_path <- args$plot_path
overwrite <- args$overwrite


# DEBUG ONLY
# mat_file <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/new/UCLA_time_series_four_groups.mat"
# label_metadata <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/participants.csv"
# rdata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/Rdata/"
# plot_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/plots/"
# overwrite <- T

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Load data from matlab into R
load_mat_data(mat_file = mat_file,
              label_metadata = label_metadata,
              rdata_path = rdata_path,
              overwrite = F)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Run catch22 on time-series data for each noise processing method
for (noise_proc in c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")) {
  noise_label <- gsub("\\+", "_", noise_proc)
  TS_df <- readRDS(paste0(rdata_path, sprintf("UCLA_%s.Rds", noise_label)))
  if (!file.exists(paste0(rdata_path, sprintf("UCLA_%s_catch22.Rds", 
                                              noise_label)))) {
    cat("\nNow running catch22 for UCLA", noise_proc, "data.\n")
    TS_catch22 <- catch22_all_regions(TS_df=TS_df)
    saveRDS(TS_catch22, file=paste0(rdata_path, sprintf("UCLA_%s_catch22.Rds", 
                                                        noise_label)))
  }
  # clean up memory
  gc()
  remove(TS_df)
}

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Run univariate classification on each brain region separately
region_wise_univ_class_res_list <- list()
for (region in unique(UCLA_AROMA_2P_catch22$Brain_Region)) {
  region_res <- region_by_region_analysis(ROI=region,
                                          test_method="t-test",
                                          feature_matrix=UCLA_AROMA_2P_catch22,
                                          display_figures = F,
                                          return_restable = T)
  region_res$Brain_Region <- region
  region_wise_univ_class_res_list[[region]] <- region_res
}
region_wise_univ_class_res <- do.call(plyr::rbind.fill, region_wise_univ_class_res_list)

# TO-DO: histogram of accuracies across all regions a la Figure 1C in old readme
theme_set(cowplot::theme_cowplot())
region_wise_univ_class_res %>%
  # filter(feature %in% unique(region_wise_univ_class_res$feature)[1:2]) %>%
  mutate(feature = str_replace_all(feature, "_", " ")) %>%
  mutate(feature = str_replace_all(feature, "catch22 ", "")) %>%
  ggplot(data=., mapping=aes(x=statistic_value)) +
  geom_histogram(fill="lightsteelblue") +
  geom_vline(xintercept=0, linetype=2) +
  xlab("T statistic") +
  ylab("Number of ROIs") +
  facet_wrap(feature ~ ., scales="free", nrow=4,
             labeller = labeller(feature = label_wrap_gen(26)))
ggsave(paste0(plot_path, "catch22_T_score_histograms.png"),
       width=15, height=10, units="in", dpi=300)