# Load needed libraries
library(tidyverse)

# Define paths specific to the UCLA CNP and ABIDE ASD datasets
univariate_feature_set <- "catch22"
pairwise_feature_set <- "pyspi14"
github_dir <- "~/github/fMRI_FeaturesDisorders/"
data_path <- "~/data/UCLA_CNP_ABIDE_ASD/"
rdata_path <- paste0(data_path, "processed_data/Rdata/")
UCLA_CNP_noise_proc <- "AROMA+2P+GMR"
ABIDE_ASD_noise_proc <- "FC1000"

############################## Metadata ########################################

# Load UCLA CNP sample metadata
UCLA_CNP_metadata <- readRDS(paste0(data_path, "study_metadata/UCLA_CNP_sample_metadata.Rds")) %>%
  dplyr::rename("Sex" = "gender", "Age" = "age") %>%
  mutate(Study = "UCLA_CNP") %>%
  filter(!is.na(Diagnosis)) %>%
  mutate(Diagnosis = ifelse(Diagnosis == "ADHD", "ADHD", str_to_title(Diagnosis)))

# Load ABIDE ASD sample metadata
ABIDE_ASD_metadata <- readRDS(paste0(data_path, "study_metadata/ABIDE_ASD_sample_metadata.Rds")) %>%
  dplyr::rename("Sex" = "sex", "Age" = "age", "Site" = "site") %>%
  dplyr::select(-asd) %>%
  mutate(Study = "ABIDE_ASD") %>%
  filter(!is.na(Diagnosis))

UCLA_CNP_ABIDE_ASD_metadata <- plyr::rbind.fill(UCLA_CNP_metadata,
                                                ABIDE_ASD_metadata)

# Save to a joint Rds file
saveRDS(UCLA_CNP_ABIDE_ASD_metadata, 
        file=paste0(data_path,
                    "study_metadata/UCLA_CNP_ABIDE_ASD_sample_metadata.Rds"))

############################## Movement ########################################

# Load UCLA CNP movement data
UCLA_CNP_movement <- read.table(paste0(data_path, "movement_data/UCLA_CNP/UCLA_CNP_mFD.txt"),
                                sep=",")
colnames(UCLA_CNP_movement) <- c("Sample_ID", "Jenkinson", "Power", "VanDijk")
UCLA_CNP_movement$Study <- "UCLA_CNP"

# Load ABIDE ASD movement data
ABIDE_ASD_movement <- read.table(paste0(data_path, "movement_data/ABIDE_ASD/ABIDE_ASD_mFD.txt"),
                                 sep=",", colClasses = c("V1" = "character"))
colnames(ABIDE_ASD_movement) <- c("Sample_ID", "Jenkinson", "Power", "VanDijk")
ABIDE_ASD_movement$Study <- "ABIDE_ASD"

# Merge the data and joint with metadata
UCLA_CNP_ABIDE_ASD_movement <- plyr::rbind.fill(UCLA_CNP_movement,
                                                ABIDE_ASD_movement) %>%
  left_join(., UCLA_CNP_ABIDE_ASD_metadata)

# Save to an Rds file
saveRDS(UCLA_CNP_ABIDE_ASD_movement, 
        file=paste0(data_path, 
                    "movement_data/UCLA_CNP_ABIDE_ASD_movement_FD.Rds"))

############################# Univariate #######################################

# Load UCLA CNP filtered and z-scored catch22 results
UCLA_CNP_catch22 <- readRDS(paste0(rdata_path,
                                   "UCLA_CNP_catch22_filtered_zscored.Rds")) %>%
  filter(Noise_Proc == UCLA_CNP_noise_proc) %>%
  mutate(Study = "UCLA_CNP",
         feature_set = univariate_feature_set)

# Load UCLA CNP filtered catch2 results
UCLA_CNP_catch2 <- readRDS(paste0(rdata_path, "UCLA_CNP_catch2_filtered.Rds")) %>%
  filter(Noise_Proc == UCLA_CNP_noise_proc) %>%
  filter(Sample_ID %in% unique(UCLA_CNP_catch22$Sample_ID)) %>%
  mutate(Study = "UCLA_CNP",
         feature_set = "catch2")

# Load ABIDE ASD filtered and z-scored catch22 results
ABIDE_ASD_catch22 <- readRDS(paste0(rdata_path,
                                   "ABIDE_ASD_catch22_filtered_zscored.Rds")) %>%
  filter(Noise_Proc == ABIDE_ASD_noise_proc) %>%
  mutate(Study = "ABIDE_ASD",
         feature_set = univariate_feature_set)

# Load ABIDE ASD filtered catch2 results
ABIDE_ASD_catch2 <- readRDS(paste0(rdata_path, "ABIDE_ASD_catch2_filtered.Rds")) %>%
  filter(Noise_Proc == ABIDE_ASD_noise_proc) %>%
  filter(Sample_ID %in% unique(ABIDE_ASD_catch22$Sample_ID)) %>%
  mutate(Study = "ABIDE_ASD",
         feature_set = "catch2")

# Merge the univariate results
UCLA_CNP_ABIDE_ASD_univariate_features <- do.call(plyr::rbind.fill,
                                                  list(UCLA_CNP_catch22,
                                                       UCLA_CNP_catch2,
                                                       ABIDE_ASD_catch22,
                                                       ABIDE_ASD_catch2)) %>%
  dplyr::select(-method)

# Save merged univariate features to an Rds file
saveRDS(UCLA_CNP_ABIDE_ASD_univariate_features, 
        file=paste0(rdata_path, "UCLA_CNP_ABIDE_ASD_", univariate_feature_set,
                    "_and_catch2_filtered_zscored.Rds"))

############################ Pairwise ##########################################

# Load UCLA CNP filtered and z-scored pyspi14 results (when they're done)


# Load ABIDE ASD filtered and z-scored pyspi14 results (when they're done)


################# Merge subjects with univariate + pairwise data ###############
univariate_samples <- UCLA_CNP_ABIDE_ASD_univariate_features %>%
  distinct(Sample_ID) %>%
  pull(Sample_ID)
  
saveRDS(univariate_samples,
        paste0(rdata_path, "UCLA_CNP_ABIDE_ASD_samples_with_univariate_",
               univariate_feature_set,
               "_and_pairwise_", pairwise_feature_set,
               "_filtered.Rds"))

################# Organise the analysis contrasts per study ###############

study_comparisons_to_control <- UCLA_CNP_ABIDE_ASD_metadata %>%
  filter(Sample_ID %in% univariate_samples) %>%
  distinct(Study, Diagnosis) %>%
  rowwise() %>%
  filter(Diagnosis != "Control") %>%
  dplyr::rename("Group_to_Compare" = "Diagnosis")

saveRDS(study_comparisons_to_control,
        paste0(rdata_path, "UCLA_CNP_ABIDE_ASD_comparisons_to_control.Rds"))