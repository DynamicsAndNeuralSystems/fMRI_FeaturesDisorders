################################################################################
# Load libraries
################################################################################

library(tidyverse)
library(knitr)
library(kableExtra)

################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
data_path <- "~/data/UCLA_CNP_ABIDE_ASD/"
rdata_path <- paste0(data_path, "processed_data/Rdata/")

study_metadata <- readRDS(paste0(data_path, "study_metadata/UCLA_CNP_ABIDE_ASD_sample_metadata.Rds"))

study_metadata %>%
  group_by(Study, Diagnosis) %>%
  summarise(N = n(),
            N_Female = sum(Sex=="F"),
            Percent_Female = round(100*sum(Sex=="F")/N, 1),
            Age_Mean = round(mean(Age), 1),
            Age_SD = round(sd(Age), 1)) %>%
  rowwise() %>%
  mutate(`% Female (N)` = paste0(Percent_Female, " (", N_Female, ")"),
         `Age; Mean (SD)` = paste0(Age_Mean, " (", Age_SD, ")"),
         .keep = "unused") %>%
  kable() %>%
  kable_styling(full_width=F)
  