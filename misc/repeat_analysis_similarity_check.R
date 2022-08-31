library(tidyverse)

data_path <- "/headnode1/abry4213/data/UCLA_Schizophrenia/"

# Load in data for each of 5 catch22 runs after filtering and z-scoring
catch22_list <- list()
for (run_number in 1:5){
  run_number_catch22 <- readRDS(sprintf("%s/processed_data_run%s/Rdata/UCLA_Schizophrenia_catch22_filtered_zscored.Rds",
                                        data_path, run_number))
  catch22_list[[run_number]] <- run_number_catch22
}

# Check whether all five dataframes in terms are equal
all(apply(combn(length(catch22_list), 2), 2, function(x)
  all.equal(catch22_list[[x[1]]], catch22_list[[x[2]]])))
# TRUE woo

