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
                                unique_columns = c("Sample_ID", "Brain_Region", "Noise_Proc"),
                                overwrite = F) {
  
  # Save resulting time-series features to an .Rds file
  if (!file.exists(paste0(rdata_path, sprintf("%s_catch22.Rds", 
                                              input_dataset_name))) | overwrite) {
    
    # Load data
    TS_df <- readRDS(TS_data_file)
    
    # Merge columns into unique ID
    TS_df <- TS_df %>%
      tidyr::unite("Unique_ID", unique_columns, sep="__")
    
    # Run catch22 using theft
    TS_catch22 <- theft::calculate_features(data = TS_df, 
                                            id_var = "Unique_ID", 
                                            time_var = "timepoint", 
                                            values_var = "value", 
                                            feature_set = "catch22",
                                            catch24 = F) %>%
      tidyr::separate("id", c(output_column_names), sep="__")
    
    cat("\nNow saving resulting theft features to an .Rds file: \n", 
        paste0(rdata_path, sprintf("%s_catch22.Rds", 
                                   input_dataset_name)))
    saveRDS(TS_catch22, file=paste0(rdata_path, sprintf("%s_catch22.Rds", 
                                                        input_dataset_name)))
  }
}
