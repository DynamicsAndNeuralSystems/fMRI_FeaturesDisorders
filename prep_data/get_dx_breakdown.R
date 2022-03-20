#------------------------------------
# This script sets out to produce a
# function for reading in matlab time
# series files into R
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 15 March 2022
#--------------------------------------

library(tidyverse)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Get breakdown of control vs schizophrenia subjects in cohort

get_dx_breakdown <- function(subject_csv) {
  dx_data <- read.csv(subject_csv) %>%
    filter(diagnosis %in% c("Control", "Schz"))
  num_total <- nrow(dx_data)
  num_scz <- nrow(filter(dx_data, diagnosis=="Schz"))
  num_ctrl <- nrow(filter(dx_data, diagnosis=="Control"))
  scz2ctrl <- num_scz/num_ctrl
  
  res_df <- data.frame(var1 = num_ctrl,
                       var2 = num_scz,
                       var3 = num_total,
                       var4 = scz2ctrl)
  colnames(res_df) <- c("Control", "Schizophrenia", "Total", "SCZ2Ctrl")
  print(res_df)
}