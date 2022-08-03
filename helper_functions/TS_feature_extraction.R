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

#-------------------------------------------------------------------------------
# Helper function to run catch22 with theft for all subject/brain region combos
#-------------------------------------------------------------------------------

catch22_all_samples <- function(TS_data_file,
                                rdata_path,
                                input_dataset_name = "UCLA",
                                output_column_names = c("Output"),
                                overwrite = F) {
  # Load data
  TS_df <- readRDS(TS_data_file)
  
  # Run catch22 using theft
  TS_catch22 <- theft::calculate_features(data = TS_df, 
                                          id_var = "Unique_ID", 
                                          time_var = "timepoint", 
                                          values_var = "value", 
                                          feature_set = "catch22",
                                          catch24 = F) %>%
    tidyr::separate("id", c(output_column_names), sep="__")
  
  # Save resulting time-series features to an .Rds file
  if (!file.exists(paste0(rdata_path, sprintf("%s_catch22.Rds", 
                                              input_dataset_name))) | overwrite) {
    cat("\nNow saving resulting theft features to an .Rds file: \n", 
        paste0(rdata_path, sprintf("%s_catch22.Rds", 
                                   input_dataset_name)))
    saveRDS(TS_catch22, file=paste0(rdata_path, sprintf("%s_catch22.Rds", 
                                                        input_dataset_name)))
  }
}
