#------------------------------------
# This script sets out to produce a
# function for reading in matlab time
# series files into R
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 1 August 2022
#--------------------------------------

require(theft) 
require(tidyr)
require(feather)

#-------------------------------------------------------------------------------
# Helper function to run catch24 with theft for all subject/brain region combos
#-------------------------------------------------------------------------------

catch24_all_samples <- function(full_TS_data,
                                data_path,
                                dataset_ID = "UCLA_CNP",
                                noise_proc = "AROMA+2P+GMR",
                                output_column_names = c("Output"),
                                unique_columns = c("Sample_ID", "Brain_Region", "Noise_Proc")) {
  
  noise_label = gsub("\\+", "_", noise_proc)
  
  # Merge columns into unique ID
  full_TS_data <- full_TS_data %>%
    tidyr::unite("Unique_ID", unique_columns, sep="__")
  
  # Run catch24 using theft
  TS_catch24 <- theft::calculate_features(data = full_TS_data, 
                                          id_var = "Unique_ID", 
                                          time_var = "timepoint", 
                                          values_var = "values", 
                                          feature_set = "catch22",
                                          catch24 = TRUE)[[1]] %>%
    tidyr::separate("id", c(output_column_names), sep="__")
  
  cat("\nNow saving resulting theft features to a .feather file: \n", 
      glue("{data_path}/raw_data/{dataset_ID}_{noise_label}_catch24.feather"))

    # Save full catch24 feature set
    feather::write_feather(TS_catch24, glue("{data_path}/raw_data/{dataset_ID}_{noise_label}_catch24.feather"))
    
}
