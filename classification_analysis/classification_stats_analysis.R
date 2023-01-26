#-------------------------------------------------------------------------------
# Define paths
#-------------------------------------------------------------------------------

python_to_use <- "~/.conda/envs/pyspi/bin/python3"
# python_to_use <- "/Users/abry4213/opt/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
github_dir <- "~/github/"
data_path <- "~/data/TS_feature_manuscript"
study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Noise_Proc = c(rep("AROMA+2P+GMR",3), "FC1000"),
                             Comparison_Group = c("Schizophrenia", "ADHD", "Bipolar", "ASD"))
univariate_feature_sets <- c("catch22", "catch2", "catch24")

list.append <- function (.data, ...) 
{
  if (is.list(.data)) {
    c(.data, list(...))
  }
  else {
    c(.data, ..., recursive = FALSE)
  }
}

reticulate::use_python(python_to_use)
library(feather)
library(tidyverse)
library(reticulate)
library(glue)

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

# Load p-value results and null distributions
univariate_p_values <- pyarrow_feather$read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_empirical_p_values.feather"))
pairwise_p_values <- pyarrow_feather$read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_empirical_p_values.feather"))
combo_univariate_pairwise_p_values <- pyarrow_feather$read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_combo_univariate_pairwise_empirical_p_values.feather"))

univariate_null_distribution <- pyarrow_feather$read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_univariate_null_balanced_accuracy_distributions.feather"))
pairwise_null_distribution <- pyarrow_feather$read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_pairwise_null_balanced_accuracy_distributions.feather"))
combo_univariate_pairwise_null_distribution <- pyarrow_feather$read_feather(glue("{output_data_path}/UCLA_CNP_ABIDE_ASD_combo_univariate_pairwise_null_balanced_accuracy_distributions.feather"))

