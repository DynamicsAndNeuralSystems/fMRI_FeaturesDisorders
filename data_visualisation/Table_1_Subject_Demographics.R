################################################################################
# Load libraries
################################################################################

python_to_use <- "/Users/abry4213/opt/anaconda3/envs/pyspi/bin/python3"
reticulate::use_python(python_to_use)
library(reticulate)

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

library(tidyverse)
library(knitr)
library(kableExtra)
library(glue)
library(rstatix)

################################################################################
# Define study/data paths
################################################################################

# Filter to just those samples used in each study
UCLA_CNP_subjects_to_keep <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_filtered_sample_info_AROMA_2P_GMR_catch24_pyspi14.feather")
ABIDE_ASD_subjects_to_keep <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/processed_data/ABIDE_ASD_filtered_sample_info_FC1000_catch24_pyspi14.feather")

UCLA_CNP_study_metadata <- pyarrow_feather$read_feather("~/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather") %>%
  semi_join(., UCLA_CNP_subjects_to_keep) %>%
  mutate(Age = as.numeric(Age))
ABIDE_ASD_study_metadata <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_sample_metadata.feather") %>%
  semi_join(., ABIDE_ASD_subjects_to_keep) %>%
  mutate(Age = as.numeric(Age))

# Generate Table 1
plyr::rbind.fill(UCLA_CNP_study_metadata, ABIDE_ASD_study_metadata) %>%
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
  
# Chi-square analysis for sex
# UCLA CNP Control vs. Schizophrenia
UCLA_CNP_study_metadata %>%
  filter(Diagnosis %in% c("Control", "Schizophrenia")) %>%
  count(Diagnosis, Sex) %>%
  spread(Sex, n) %>%
  select(-Diagnosis) %>%
  chisq.test()

# UCLA CNP Control vs. Bipolar
UCLA_CNP_study_metadata %>%
  filter(Diagnosis %in% c("Control", "Bipolar")) %>%
  count(Diagnosis, Sex) %>%
  spread(Sex, n) %>%
  select(-Diagnosis) %>%
  chisq.test()

# UCLA CNP Control vs. ADHD
UCLA_CNP_study_metadata %>%
  filter(Diagnosis %in% c("Control", "ADHD")) %>%
  count(Diagnosis, Sex) %>%
  spread(Sex, n) %>%
  select(-Diagnosis) %>%
  chisq.test()

# ABIDE Control vs ASD
ABIDE_ASD_study_metadata %>%
  filter(Diagnosis %in% c("Control", "ASD")) %>%
  count(Diagnosis, Sex) %>%
  spread(Sex, n) %>%
  select(-Diagnosis) %>%
  chisq.test()

# Wilcoxon rank-sum test for age
# UCLA CNP Control vs. Schizophrenia
UCLA_CNP_study_metadata %>%
  filter(Diagnosis %in% c("Control", "Schizophrenia")) %>%
  wilcox_test(Age ~ Diagnosis)

# UCLA CNP Control vs. Bipolar
UCLA_CNP_study_metadata %>%
  filter(Diagnosis %in% c("Control", "Bipolar"))%>%
  wilcox_test(Age ~ Diagnosis)

# UCLA CNP Control vs. ADHD
UCLA_CNP_study_metadata %>%
  filter(Diagnosis %in% c("Control", "ADHD")) %>%
  wilcox_test(Age ~ Diagnosis)

# ABIDE Control vs ASD
ABIDE_ASD_study_metadata %>%
  filter(Diagnosis %in% c("Control", "ASD")) %>%
  wilcox_test(Age ~ Diagnosis)