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
# Helper functions to run catch22 with theft for all subject/brain region combos
#-------------------------------------------------------------------------------

catch22_by_noise_proc <- function(noise_proc,
                                  rdata_path,
                                  TS_df,
                                  input_dataset_name = "UCLA") {
  
  # Clean up noise processing label for file naming
  noise_label <- gsub("\\+", "_", noise_proc)
  
  # Filter time-series data to noise processing type
  TS_df <- dplyr::filter(TS_df, Noise_Proc == noise_proc)
  
  # Create unique ID variable that combines subject ID with brain region
  TS_df$unique_ID <- paste(TS_df$Subject_ID, TS_df$Brain_Region, sep="_")
  
  # Run catch22 using theft
  TS_catch22 <- theft::calculate_features(data = TS_df, 
                                          id_var = "unique_ID", 
                                          time_var = "timepoint", 
                                          values_var = "value", 
                                          group_var = "diagnosis", 
                                          feature_set = "catch22",
                                          catch24 = F) %>%
    tidyr::separate(id, into=c("Subject_ID", "Brain_Region"), sep="_") %>%
    dplyr::mutate(Noise_Proc = noise_proc)
  
  return(TS_catch22)
}

catch22_all_regions <- function(TS_data_file,
                                rdata_path,
                                input_dataset_name = "UCLA",
                                noise_procs = c("AROMA+2P",
                                                "AROMA+2P+GMR", 
                                                "AROMA+2P+DiCER")) {
  # Load data
  TS_df <- readRDS(TS_data_file)
  
  # Apply catch22_by_noise_proc to each noise-processing method
  catch22_res <- noise_procs %>%
    purrr::map_df(~ catch22_by_noise_proc(noise_proc = .x,
                                          TS_df = TS_df,
                                          rdata_path = rdata_path,
                                          input_dataset_name = input_dataset_name))
  
  # Save resulting time-series features to an .Rds file
  saveRDS(catch22_res, file=paste0(rdata_path, sprintf("%s_catch22.Rds", 
                                                      input_dataset_name)))
}
