#------------------------------------
# This script sets out to produce a
# function for reading in matlab time
# series files into R
#------------------------------------

#--------------------------------------
# Author: Trent Henderson, 9 March 2021
# Updated: Annie Bryant, 15 March 2022
#--------------------------------------

library(tidyverse)
library(theft)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# catch22 for all regions

catch22_all_regions <- function(TS_df) {
  # Create a new ID that contains both subject ID and brain region
  TS_df$unique_ID <- paste(TS_df$Subject_ID, TS_df$Brain_Region, sep="_")
  
  # Run catch22 using theft
  feature_matrix <- calculate_features(data = TS_df, 
                                       id_var = "unique_ID", 
                                       time_var = "timepoint", 
                                       values_var = "value", 
                                       group_var = "diagnosis", 
                                       feature_set = "catch22",
                                       catch24 = F) %>%
    tidyr::separate(id, into=c("Subject_ID", "Brain_Region"), sep="_")
  
  return(feature_matrix)
}
