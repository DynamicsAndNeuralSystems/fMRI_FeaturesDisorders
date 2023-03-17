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
# Helper function to run catch22 with theft for all subject/brain region combos
#-------------------------------------------------------------------------------

catch22_all_samples <- function(full_TS_data,
                                rdata_path,
                                dataset_ID = "UCLA_CNP",
                                noise_proc = "AROMA+2P+GMR",
                                output_column_names = c("Output"),
                                unique_columns = c("Sample_ID", "Brain_Region", "Noise_Proc"),
                                add_mean_SD = FALSE,
                                overwrite = F) {

  noise_label = gsub("\\+", "_", noise_proc)
  
  # Save resulting time-series features to an .Rds file
  if (!file.exists(paste0(rdata_path, sprintf("%s_catch22.feather", 
                                              dataset_ID))) | overwrite) {
    
    # Merge columns into unique ID
    full_TS_data <- full_TS_data %>%
      tidyr::unite("Unique_ID", unique_columns, sep="__")
    
    # Run catch22 using theft
    TS_catch22 <- theft::calculate_features(data = full_TS_data, 
                                            id_var = "Unique_ID", 
                                            time_var = "timepoint", 
                                            values_var = "values", 
                                            feature_set = "catch22",
                                            catch24 = add_mean_SD)[[1]] %>%
      tidyr::separate("id", c(output_column_names), sep="__")
    
    cat("\nNow saving resulting theft features to a .feather file: \n", 
        paste0(rdata_path, sprintf("%s_%s_catch22.feather", 
                                   dataset_ID, noise_label)))
    if (add_mean_SD) {
      # Save full catch24 feature set
      TS_catch24 <- TS_catch22
      feather::write_feather(TS_catch24, paste0(rdata_path, sprintf("%s_%s_catch24.feather",
                                                          dataset_ID, noise_label)))
      
      # Also save just catch22
      TS_catch22 <- TS_catch24 %>%
        filter(!(names %in% c("DN_Mean", "DN_Spread_Std")))
      
      feather::write_feather(TS_catch22, paste0(rdata_path, sprintf("%s_%s_catch22.feather", 
                                                          dataset_ID, noise_label)))
      
      # Also save just catch2
      TS_catch2 <- TS_catch24 %>%
        filter(names %in% c("DN_Mean", "DN_Spread_Std"))
      feather::write_feather(TS_catch2, paste0(rdata_path, sprintf("%s_%s_catch2.feather", 
                                                         dataset_ID, noise_label)))
    } else {
      feather::write_feather(TS_catch22, paste0(rdata_path, sprintf("%s_%s_catch22.feather", 
                                                          dataset_ID, noise_label)))
    }
  }
}
