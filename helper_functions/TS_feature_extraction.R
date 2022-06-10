#------------------------------------
# This script sets out to produce a
# function for reading in matlab time
# series files into R
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 8 June 2022
#--------------------------------------

require(theft) 
require(tidyr)

#-------------------------------------------------------------------------------
# Helper function to run catch22 with theft for all subject/brain region combos
#-------------------------------------------------------------------------------

catch22_all_regions <- function(rdata_path,
                                input_dataset_name = "UCLA",
                                noise_procs = c("AROMA+2P",
                                                "AROMA+2P+GMR", 
                                                "AROMA+2P+DiCER")) {
  for (noise_proc in c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")) {
    
    # Clean up noise processing label for file naming
    noise_label <- gsub("\\+", "_", noise_proc)
    
    # Create Rds file for given noise processing method if it doesn't already exist
    if (!file.exists(paste0(rdata_path, sprintf("%s_%s_catch22.Rds", 
                                                input_dataset_name, noise_label)))) {
      
      # Read in time-series data
      TS_df <- readRDS(paste0(rdata_path, sprintf("UCLA_%s.Rds", noise_label)))
      cat("\nNow running catch22 for", input_dataset_name, noise_proc, "data.\n")
      
      # Create unique ID that combines subject ID with brain region
      TS_df$unique_ID <- paste(TS_df$Subject_ID, TS_df$Brain_Region, sep="_")
      
      # Run catch22 using theft
      TS_catch22 <- theft::calculate_features(data = TS_df, 
                                              id_var = "unique_ID", 
                                              time_var = "timepoint", 
                                              values_var = "value", 
                                              group_var = "diagnosis", 
                                              feature_set = "catch22",
                                              catch24 = F) %>%
        tidyr::separate(id, into=c("Subject_ID", "Brain_Region"), sep="_")
      
      # Save resulting time-series features to an .Rds file
      saveRDS(TS_catch22, file=paste0(rdata_path, sprintf("%s_%s_catch22.Rds", 
                                                          input_dataset_name,
                                                          noise_label)))
      
    }
  }
}
