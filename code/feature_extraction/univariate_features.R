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
                                output_column_names = c("Output"),
                                unique_columns = c("Sample_ID", "Brain_Region", "Noise_Proc")) {
  
  
  # Merge columns into unique ID
  full_TS_data <- full_TS_data %>%
    tidyr::unite("Unique_ID", unique_columns, sep="__")
  
  # Compute the set of 24 time-series features using theft
  TS_catch24 <- theft::calculate_features(data = full_TS_data, 
                                          id_var = "Unique_ID", 
                                          time_var = "timepoint", 
                                          values_var = "values", 
                                          feature_set = "catch22",
                                          catch24 = TRUE)[[1]] %>%
    tidyr::separate("id", c(output_column_names), sep="__")

  # Return the resulting set of 24 features computed per brain region
  return(TS_catch24)
    
}
