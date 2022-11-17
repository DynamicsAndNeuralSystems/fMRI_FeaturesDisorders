################################################################################
# Load libraries
################################################################################

library(tidyverse)
library(icesTAF)
library(cowplot)
theme_set(theme_cowplot())

################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/")
icesTAF::mkdir(plot_path)

source(paste0(github_dir, "helper_functions/data_prep_and_QC/QC_functions_univariate.R"))
source(paste0(github_dir, "helper_functions/Visualization.R"))

UCLA_data_path <- "~/data/UCLA_Schizophrenia/"
UCLA_rdata_path <- paste0(UCLA_data_path, "processed_data/Rdata/")
ABIDE_data_path <- "~/data/ABIDE_ASD/"
ABIDE_rdata_path <- paste0(ABIDE_data_path, "processed_data/Rdata/")


ABIDE_brain_region_info <- read.csv(paste0(ABIDE_data_path, 
                                           "Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv"))

icesTAF::mkdir(paste0(plot_path, "Figure2/"))
