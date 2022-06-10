# TODO: convert to argparse
python_to_use <- "/home/osboxes/anaconda3/envs/pyspi/bin/python3"
reticulate::use_python(python_to_use)
set.seed(127)

github_dir <- "/media/sf_Shared_Folder/github/fMRI_FeaturesDisorders/"

study <- "/media/sf_Shared_Folder/PhD_work/"
data_path <- paste0(study, "data/scz/UCLA/")
pydata_path <- paste0(study, "data/scz/UCLA/pydata/")
output_data_path <- paste0(study, "data/scz/UCLA/pydata/R_files/")

# load libraries
library(theft)
library(tidyverse)
library(cowplot)
library(reticulate)
source_python(paste0(github_dir, "helper_functions/pickle_reader.py"))
theme_set(theme_cowplot())

# Load subject metadata
subject_csv <- read.csv(paste0(data_path, "participants.csv"))

################################################################################
# UCLA data prep -- LOCAL ON UBUNTU
################################################################################

### Source helper function
# TODO: This python script needs to have argparse added
system(sprintf("python3 %s/helper_functions/split_MTS_into_npy.py",
               github_dir))


################################################################################
# Run pySPI with reduced SPI set -- ON PHYSICS CLUSTER
################################################################################

### See github_dir/pyspi_files/call_pyspi_on_cluster.sh

################################################################################
# Read pyspi calc.pkl files into R for each subject -- LOCAL ON UBUNTU
################################################################################
noise_procs = c("AROMA+2P",
                "AROMA+2P+GMR",
                "AROMA+2P+DiCER")

source(paste0(github_dir, "helper_functions/pyspi_functions.R"))
read_pyspi_pkl_into_RDS(data_path = pydata_path,
                        noise_procs = noise_procs)

# TODO: quality control for pyspi data