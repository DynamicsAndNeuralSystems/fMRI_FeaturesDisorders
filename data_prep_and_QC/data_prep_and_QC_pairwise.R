# TODO: convert to argparse
python_to_use <- "/home/osboxes/anaconda3/envs/pyspi/bin/python3"

# github_dir <- "/media/sf_Shared_Folder/github/fMRI_FeaturesDisorders/"
# study <- "/media/sf_Shared_Folder/PhD_work/"
github_dir <- "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/"
study <- "D:/Virtual_Machines/Shared_Folder/PhD_work/"
data_path <- paste0(study, "data/scz/UCLA/")
rdata_path <- paste0(data_path, "Rdata/")
pydata_path <- paste0(data_path, "pydata/")
output_data_path <- paste0(study, "data/scz/UCLA/pydata/R_files/")

# load libraries
reticulate::use_python(python_to_use)
library(theft)
library(tidyverse)
library(cowplot)
library(reticulate)
theme_set(theme_cowplot())
source_python(paste0(github_dir, "helper_functions/pickle_reader.py"))

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
feature_set = "pyspi_19"

################################################################################
# Read pyspi calc.pkl files into R for each subject -- LOCAL ON UBUNTU
################################################################################
noise_procs = c("AROMA+2P",
                "AROMA+2P+GMR",
                "AROMA+2P+DiCER")

source(paste0(github_dir, "helper_functions/pyspi_functions.R"))
read_pyspi_pkl_into_RDS(data_path = pydata_path,
                        subject_csv_file = paste0(data_path, "participants.csv"),
                        noise_procs = noise_procs)

################################################################################
# Merge pyspi res for each subject -- LOCAL ON UBUNTU
################################################################################
input_dataset_name = "UCLA"
ROI_index_file <- paste0(study, "data/scz/UCLA/pydata/ROI_info.csv")
merge_pyspi_res_for_study(data_path = pydata_path,
                          input_dataset_name = input_dataset_name,
                          ROI_index_file = ROI_index_file,
                          noise_procs = noise_procs)

################################################################################
# QC for pyspi res
################################################################################

aroma_2P_GMR_pyspi <- readRDS(paste0(pydata_path, "UCLA_all_subject_pyspi_AROMA_2P_GMR.Rds"))

### Number of NA Values by ROI Pair
aroma_2P_GMR_pyspi %>%
  tidyr::unite("ROI_pair", c(brain_region_1, brain_region_2), sep="_") %>%
  group_by(ROI_pair, SPI) %>%
  summarise(num_subj_NA = sum(is.na(value))) %>%
  ggplot(data=., mapping = aes(y=SPI, x=ROI_pair, fill=num_subj_NA)) +
  ggtitle("Number of NA Values per\nROI Pair + SPI") +
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  labs(fill = "# Subjects") +
  xlab("ROI Pair") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust=0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
  

## Number of NAs by Subject
aroma_2P_GMR_pyspi %>%
  tidyr::unite("ROI_pair", c(brain_region_1, brain_region_2), sep="_") %>%
  group_by(Subject_ID, SPI) %>%
  summarise(num_ROI_NA = sum(is.na(value))) %>%
  ggplot(data=., mapping = aes(y=SPI, x=Subject_ID, fill=num_ROI_NA)) +
  ggtitle("Number of NA Values\nper Subject + SPI") +
  xlab("Subject ID") +
  labs(fill = "# ROI Pairs with NA") +
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust=0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

## Remove subject with all NA values
pyspi_all_NA_subjects <- aroma_2P_GMR_pyspi %>%
  tidyr::unite("ROI_pair", c(brain_region_1, brain_region_2), sep="_") %>%
  group_by(Subject_ID) %>%
  filter(all(is.na(value))) %>%
  ungroup() %>%
  dplyr::distinct(Subject_ID) %>%
  pull(Subject_ID)
pyspi_all_NA_subjects

## Load original TS data and plot raw time-series for NA subject
original_TS_data <- readRDS(paste0(data_path, "Rdata/UCLA_AROMA_2P_GMR.Rds"))

original_TS_data %>%
  filter(Subject_ID %in% pyspi_all_NA_subjects) %>%
  ggplot(data=., mapping=aes(x=timepoint, y=value, 
                             color=Brain_Region, group=Brain_Region)) +
  geom_line() +
  facet_wrap(Subject_ID ~ ., scales="free") +
  theme(legend.position="none")

### remove NA subjects
aroma_2P_GMR_pyspi_filtered <- aroma_2P_GMR_pyspi %>%
  filter(!(Subject_ID %in% pyspi_all_NA_subjects))

saveRDS(aroma_2P_GMR_pyspi_filtered, file=paste0(pydata_path, "UCLA_all_subject_pyspi_AROMA_2P_GMR_filtered.Rds"))

# Save subject info
filtered_subject_info <- aroma_2P_GMR_pyspi_filtered %>%
  distinct(Subject_ID, group)
saveRDS(filtered_subject_info, file=paste0(rdata_path, sprintf("Filtered_subject_info_%s.Rds",
                                                               feature_set)))

################################################################################
# Z-score normalisation
################################################################################
### Source helper function
source(paste0(github_dir, "helper_functions/QC_functions.R"))

# TODO: z-score normalise filtered pyspi data
aroma_2P_GMR_pyspi_filtered_z <- aroma_2P_GMR_pyspi_filtered %>%
  group_by(SPI) %>%
  mutate(value_z = (value - mean(value, na.rm=T))/sd(value, na.rm=T)) %>%
  dplyr::select(-value) %>%
  dplyr::rename("value" = "value_z")

saveRDS(aroma_2P_GMR_pyspi_filtered_z, file = paste0(pydata_path, "UCLA_all_subject_pyspi_AROMA_2P_GMR_filtered_zscored.Rds"))
