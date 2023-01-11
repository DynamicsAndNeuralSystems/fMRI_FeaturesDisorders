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

UCLA_CNP_noise_label <- gsub("\\+", "_", UCLA_CNP_noise_proc)
ABIDE_ASD_noise_label <- gsub("\\+", "_", ABIDE_ASD_noise_proc)

############################## Metadata ########################################

if (!file.exists(paste0(data_path,
                        "study_metadata/UCLA_CNP_ABIDE_ASD_sample_metadata.Rds"))) {
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
} else {
  UCLA_CNP_ABIDE_ASD_metadata <- readRDS(paste0(data_path, 
                                                "study_metadata/UCLA_CNP_ABIDE_ASD_sample_metadata.Rds"))
}

############################## Movement ########################################

if (!file.exists(paste0(data_path, 
                        "movement_data/UCLA_CNP_ABIDE_ASD_movement_FD.Rds"))) {
  # Load UCLA CNP movement data
  UCLA_CNP_movement <- read.table(paste0(data_path, 
                                         "movement_data/UCLA_CNP/UCLA_CNP_mFD.txt"),
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
} else {
  UCLA_CNP_ABIDE_ASD_movement <- readRDS(paste0(data_path, 
                                                "movement_data/UCLA_CNP_ABIDE_ASD_movement_FD.Rds"))
}

############################# Univariate #######################################

### univariate_feature_set 
if (!file.exists(paste0(rdata_path, sprintf("UCLA_CNP_%s_ABIDE_ASD_%s_%s_filtered_zscored.Rds",
                                            UCLA_CNP_noise_label,
                                            ABIDE_ASD_noise_label,
                                            univariate_feature_set)))) {
  # Load UCLA CNP filtered and z-scored univariate_feature_set results
  UCLA_CNP_univariate_feature_set <- readRDS(paste0(rdata_path,
                                                    sprintf("UCLA_CNP_%s_%s_filtered_zscored.Rds",
                                                            UCLA_CNP_noise_label, 
                                                            univariate_feature_set))) %>%
    filter(Noise_Proc == UCLA_CNP_noise_proc) %>%
    mutate(Study = "UCLA_CNP",
           feature_set = univariate_feature_set)
  
  # Load ABIDE ASD filtered and z-scored univariate_feature_set results
  ABIDE_ASD_univariate_feature_set <- readRDS(paste0(rdata_path,
                                                     sprintf("ABIDE_ASD_%s_%s_filtered_zscored.Rds",
                                                             ABIDE_ASD_noise_label, 
                                                             univariate_feature_set))) %>%
    filter(Noise_Proc == ABIDE_ASD_noise_proc) %>%
    mutate(Study = "ABIDE_ASD",
           feature_set = univariate_feature_set)
  
  # Merge the univariate results
  UCLA_CNP_ABIDE_ASD_univariate_feature_set <- do.call(plyr::rbind.fill,
                                                       list(UCLA_CNP_univariate_feature_set,
                                                            ABIDE_ASD_univariate_feature_set)) %>%
    dplyr::select(-method)
  
  # Save merged univariate features to an Rds file
  saveRDS(UCLA_CNP_ABIDE_ASD_univariate_feature_set, 
          file=paste0(rdata_path, sprintf("UCLA_CNP_%s_ABIDE_ASD_%s_%s_filtered_zscored.Rds",
                                          UCLA_CNP_noise_label,
                                          ABIDE_ASD_noise_label,
                                          univariate_feature_set)))
} else {
  UCLA_CNP_ABIDE_ASD_univariate_feature_set <- readRDS(paste0(rdata_path, sprintf("UCLA_CNP_%s_ABIDE_ASD_%s_%s_filtered_zscored.Rds",
                                            UCLA_CNP_noise_label,
                                            ABIDE_ASD_noise_label,
                                            univariate_feature_set)))
}

### catch2 and catch24 as needed
if (univariate_feature_set == "catch22" & !(file.exists(paste0(rdata_path, 
                                                               sprintf("UCLA_CNP_%s_ABIDE_ASD_%s_catch2_filtered_zscored.Rds", 
                                                                       UCLA_CNP_noise_label, 
                                                                       ABIDE_ASD_noise_label))))) {
  # Load UCLA CNP filtered catch2 results
  UCLA_CNP_catch2 <- readRDS(paste0(rdata_path, sprintf("UCLA_CNP_%s_catch2_filtered.Rds", 
                                                        UCLA_CNP_noise_label))) %>%
    filter(Noise_Proc == UCLA_CNP_noise_proc) %>%
    filter(Sample_ID %in% unique(UCLA_CNP_univariate_feature_set$Sample_ID)) %>%
    mutate(Study = "UCLA_CNP",
           feature_set = "catch2")
  
  # Load ABIDE ASD filtered catch2 results
  ABIDE_ASD_catch2 <- readRDS(paste0(rdata_path, sprintf("ABIDE_ASD_%s_catch2_filtered.Rds", 
                                                         ABIDE_ASD_noise_label))) %>%
    filter(Noise_Proc == ABIDE_ASD_noise_proc) %>%
    filter(Sample_ID %in% unique(ABIDE_ASD_univariate_feature_set$Sample_ID)) %>%
    mutate(Study = "ABIDE_ASD",
           feature_set = "catch2")
  
  # Merge the univariate results
  UCLA_CNP_ABIDE_ASD_catch2 <- do.call(plyr::rbind.fill,
                                       list(UCLA_CNP_catch2,
                                            ABIDE_ASD_catch2)) %>%
    dplyr::select(-method)
  # Save merged univariate features to an Rds file
  saveRDS(UCLA_CNP_ABIDE_ASD_catch2, 
          file=paste0(rdata_path, sprintf("UCLA_CNP_%s_ABIDE_ASD_%s_catch2_filtered_zscored.Rds", 
                                          UCLA_CNP_noise_label, 
                                          ABIDE_ASD_noise_label)))
  
  # Load UCLA CNP filtered catch24 results
  UCLA_CNP_catch24 <- readRDS(paste0(rdata_path, sprintf("UCLA_CNP_%s_catch24_filtered.Rds", 
                                                         UCLA_CNP_noise_label))) %>%
    filter(Noise_Proc == UCLA_CNP_noise_proc) %>%
    filter(Sample_ID %in% unique(UCLA_CNP_univariate_feature_set$Sample_ID)) %>%
    mutate(Study = "UCLA_CNP",
           feature_set = "catch24")
  
  # Load ABIDE ASD filtered catch2 results
  ABIDE_ASD_catch24 <- readRDS(paste0(rdata_path, sprintf("ABIDE_ASD_%s_catch24_filtered.Rds", 
                                                          ABIDE_ASD_noise_label))) %>%
    filter(Noise_Proc == ABIDE_ASD_noise_proc) %>%
    filter(Sample_ID %in% unique(ABIDE_ASD_univariate_feature_set$Sample_ID)) %>%
    mutate(Study = "ABIDE_ASD",
           feature_set = "catch24")
  
  # Merge the univariate results
  UCLA_CNP_ABIDE_ASD_catch24 <- do.call(plyr::rbind.fill,
                                        list(UCLA_CNP_catch24,
                                             ABIDE_ASD_catch24)) %>%
    dplyr::select(-method)
  
  # Save merged univariate features to an Rds file
  saveRDS(UCLA_CNP_ABIDE_ASD_catch2, 
          file=paste0(rdata_path, sprintf("UCLA_CNP_%s_ABIDE_ASD_%s_catch24_filtered_zscored.Rds", 
                                          UCLA_CNP_noise_label, 
                                          ABIDE_ASD_noise_label)))
  
}

############################ Pairwise ##########################################

if (!file.exists(paste0(rdata_path, sprintf("UCLA_CNP_%s_ABIDE_ASD_%s_%s_filtered_zscored.Rds",
                                            UCLA_CNP_noise_label, 
                                            ABIDE_ASD_noise_label, 
                                            pairwise_feature_set)))) {
  
  # Load UCLA CNP filtered and z-scored pairwise_feature_set results
  UCLA_CNP_pairwise_feature_set <- readRDS(paste0(rdata_path,
                                                  sprintf("UCLA_CNP_%s_%s_filtered_zscored.Rds",
                                                          UCLA_CNP_noise_label, 
                                                          pairwise_feature_set))) %>%
    filter(Noise_Proc == UCLA_CNP_noise_proc) %>%
    mutate(Study = "UCLA_CNP",
           feature_set = pairwise_feature_set)
  
  # Load ABIDE ASD filtered and z-scored pairwise_feature_set results
  ABIDE_ASD_pairwise_feature_set <- readRDS(paste0(rdata_path,
                                                   sprintf("ABIDE_ASD_%s_%s_filtered_zscored.Rds",
                                                           ABIDE_ASD_noise_label, 
                                                           pairwise_feature_set))) %>%
    filter(Noise_Proc == ABIDE_ASD_noise_proc) %>%
    mutate(Study = "ABIDE_ASD",
           feature_set = pairwise_feature_set)
  
  # Merge the univariate results
  UCLA_CNP_ABIDE_ASD_pairwise_feature_set <- do.call(plyr::rbind.fill,
                                                     list(UCLA_CNP_pairwise_feature_set,
                                                          ABIDE_ASD_pairwise_feature_set))
  
  # Save merged univariate features to an Rds file
  saveRDS(UCLA_CNP_ABIDE_ASD_pairwise_feature_set, 
          file=paste0(rdata_path, sprintf("UCLA_CNP_%s_ABIDE_ASD_%s_%s_filtered_zscored.Rds",
                                          UCLA_CNP_noise_label, 
                                          ABIDE_ASD_noise_label, 
                                          pairwise_feature_set)))
} else {
  UCLA_CNP_ABIDE_ASD_pairwise_feature_set <- readRDS(paste0(rdata_path, sprintf("UCLA_CNP_%s_ABIDE_ASD_%s_%s_filtered_zscored.Rds",
                                                                                UCLA_CNP_noise_label, 
                                                                                ABIDE_ASD_noise_label, 
                                                                                pairwise_feature_set)))
}

################# Merge subjects with univariate + pairwise data ###############

# univariate_feature_set
univariate_samples <- UCLA_CNP_ABIDE_ASD_univariate_feature_set %>%
  distinct(Sample_ID) %>%
  pull(Sample_ID)

# pairwise feature set
pairwise_samples <- UCLA_CNP_ABIDE_ASD_pairwise_feature_set %>%
  distinct(Sample_ID) %>%
  pull(Sample_ID)

# find intersection between the two
common_samples <- intersect(univariate_samples, pairwise_samples)

saveRDS(common_samples, paste0(rdata_path, sprintf("UCLA_CNP_%s_ABIDE_ASD_%s_samples_with_univariate_%s_and_pairwise_%s_filtered.Rds",
                                                   UCLA_CNP_noise_label, 
                                                   ABIDE_ASD_noise_label, 
                                                   univariate_feature_set, 
                                                   pairwise_feature_set)))

# catch24
saveRDS(common_samples, paste0(rdata_path, sprintf("UCLA_CNP_%s_ABIDE_ASD_%s_samples_with_univariate_catch24_and_pairwise_%s_filtered.Rds",
                                                   UCLA_CNP_noise_label, 
                                                   ABIDE_ASD_noise_label, 
                                                   pairwise_feature_set)))

# catch2
saveRDS(common_samples, paste0(rdata_path, sprintf("UCLA_CNP_%s_ABIDE_ASD_%s_samples_with_univariate_catch2_and_pairwise_%s_filtered.Rds",
                                                   UCLA_CNP_noise_label, 
                                                   ABIDE_ASD_noise_label, 
                                                   pairwise_feature_set)))

################# Organise the analysis contrasts per study ###############

if (!file.exists(paste0(rdata_path, "UCLA_CNP_ABIDE_ASD_comparisons_to_control.Rds"))) {
  study_comparisons_to_control <- UCLA_CNP_ABIDE_ASD_metadata %>%
    filter(Sample_ID %in% univariate_samples) %>%
    distinct(Study, Diagnosis) %>%
    rowwise() %>%
    filter(Diagnosis != "Control") %>%
    dplyr::rename("Group_to_Compare" = "Diagnosis")

  saveRDS(study_comparisons_to_control,
          paste0(rdata_path, "UCLA_CNP_ABIDE_ASD_comparisons_to_control.Rds"))
}