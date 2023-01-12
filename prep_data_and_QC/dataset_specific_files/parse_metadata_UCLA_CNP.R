# Load needed libraries
library(tidyverse)

# Define paths specific to the UCLA CNP dataset
univariate_feature_set <- "catch22"
pairwise_feature_set <- "pyspi14"
github_dir <- "~/github/fMRI_FeaturesDisorders/"
data_path <- "~/data/UCLA_CNP/"
rdata_path <- paste0(data_path, "processed_data/Rdata/")
metadata_CSV <- paste0(data_path, "study_metadata/UCLA_CNP_participants.csv")
noise_proc <- "AROMA+2P+GMR"

############################## Save to Rds ########################################


if (!file.exists(paste0(data_path,
                        "study_metadata/UCLA_CNP_sample_metadata.Rds"))) {
  # Load UCLA CNP sample metadata
  UCLA_CNP_metadata <- read.csv(metadata_CSV, colClasses = "character") %>%
    dplyr::select(participant_id:gender) %>%
    dplyr::rename("Sample_ID" = "participant_id",
                  "Diagnosis" = "diagnosis",
                  "Sex" = "gender", 
                  "Age" = "age") %>%
    mutate(Study = "UCLA_CNP") %>%
    filter(!is.na(Diagnosis)) %>%
    mutate(Diagnosis = ifelse(Diagnosis == "ADHD", "ADHD", str_to_title(Diagnosis)))
  
  # Save to a joint Rds file
  saveRDS(UCLA_CNP_metadata, 
          file=paste0(data_path,
                      "study_metadata/UCLA_CNP_sample_metadata.Rds"))
  
}

