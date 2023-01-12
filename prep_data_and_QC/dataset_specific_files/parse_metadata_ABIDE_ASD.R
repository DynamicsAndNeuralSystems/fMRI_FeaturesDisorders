# Load needed libraries
library(tidyverse)

# Define paths specific to the UCLA CNP dataset
univariate_feature_set <- "catch22"
pairwise_feature_set <- "pyspi14"
github_dir <- "~/github/fMRI_FeaturesDisorders/"
data_path <- "~/data/ABIDE_ASD/"
rdata_path <- paste0(data_path, "processed_data/Rdata/")
metadata_CSV <- paste0(data_path, "study_metadata/ABIDE_ASD_participants.csv")
noise_proc <- "FC1000"

############################## Save to Rds ########################################

if (!file.exists(paste0(data_path,
                        "study_metadata/ABIDE_ASD_sample_metadata.Rds"))) {
  # Load UCLA CNP sample metadata
  ABIDE_ASD_metadata <- read.csv(metadata_CSV, colClasses = "character") %>%
    dplyr::rename("Sample_ID" = "subject_id",
                  "ASD" = "asd",
                  "Sex" = "sex", 
                  "Site" = "site") %>%
    mutate(Study = "ABIDE_ASD",
           Diagnosis = ifelse(ASD==1, "ASD", "Control"), .keep="unused") %>%
    filter(!is.na(Diagnosis))
  
  # Save to a joint Rds file
  saveRDS(ABIDE_ASD_metadata, 
          file=paste0(data_path,
                      "study_metadata/ABIDE_ASD_sample_metadata.Rds"))
  
}

