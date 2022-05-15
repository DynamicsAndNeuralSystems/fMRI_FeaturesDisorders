#------------------------------------
# This script sets out to produce a
# function for reading in matlab time
# series files into R
#------------------------------------

#--------------------------------------
# Author: Trent Henderson, 9 March 2021
# Updated: Annie Bryant, 23 March 2022
#--------------------------------------

require(plyr)
library(tidyverse)
library(R.matlab)
library(theft)
library(cowplot) 
theme_set(theme_cowplot())

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Function to load matlab .mat data for UCLA cohort
load_mat_data <- function(mat_file, subject_csv, rdata_path, overwrite=F) {
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Read in data
  cat("\nLoading in .mat file:", mat_file, "\n")
  mat_data <- readMat(mat_file)
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Load noise processing info
  noise_proc <- reshape2::melt(mat_data$noiseOptions) %>%
    dplyr::rename("noiseOptions" = "L1",
                  "noise_proc" = "value") %>%
    distinct(noise_proc, noiseOptions)
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Reshape data
  TS_data_long <- reshape2::melt(mat_data$time.series) %>%
    dplyr::rename(timepoint = Var1,
                  ROI_Index = Var2,
                  Subject_Index = Var3,
                  noiseOptions = Var4)
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Get unique IDs to join in identifiers
  ids <- reshape2::melt(mat_data$subject.list) %>% 
    group_by(value, L1) %>% 
    distinct() %>%
    dplyr::rename(Subject_ID = value,
                  Subject_Index = L1) %>%
    dplyr::select(c(Subject_ID, Subject_Index))
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Retrieve labels and clean up diagnosis names
  subject_info <- read.csv(subject_csv) %>%
    dplyr::rename(Subject_ID = 1) %>%
    distinct(Subject_ID, diagnosis, age, gender) %>%
    semi_join(., ids) %>%
    mutate(diagnosis = str_to_title(diagnosis)) %>%
    mutate(diagnosis = ifelse(diagnosis == "Adhd", "ADHD", diagnosis)) 
  
  # Save CSV file containing list of subjects with time-series data and diagnoses
  write.csv(subject_info, paste0(rdata_path, "UCLA_subjects_with_TS_data.csv"),
            row.names = F)
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Extract ROI names
  ROI_info <- reshape2::melt(mat_data$StructNames) %>%
    distinct(Var1, value) %>%
    dplyr::rename("ROI_Index"="Var1",
                  "Brain_Region"="value") %>%
    mutate(Brain_Region = gsub(" +", "", Brain_Region))
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Bring all data together
  TS_data_full <- inner_join(TS_data_long, ids, 
                             by=c("Subject_Index"="Subject_Index")) %>%
    inner_join(., subject_info, by=c("Subject_ID"="Subject_ID")) %>%
    inner_join(., noise_proc, by=c("noiseOptions"="noiseOptions")) %>%
    inner_join(., ROI_info, by=c("ROI_Index"="ROI_Index")) %>%
    dplyr::select(-noiseOptions, -Subject_Index, -ROI_Index) %>%
    filter(diagnosis %in% c("Schz", "Control"))
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Split data by noise processing
  TS_data_split <- split(TS_data_full, TS_data_full$noise_proc)
  
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Save each dataframe to an Rds object
  for (noise_proc in names(TS_data_split)) {
    tmp <- TS_data_split[[noise_proc]]
    noise_label <- gsub("\\+", "_", noise_proc)
    if (!file.exists(paste0(rdata_path, 
                            sprintf("UCLA_%s.Rds", noise_label))) | 
        overwrite) {
      cat("\nWriting UCLA", noise_proc, "data to Rds object.", "\n")
      saveRDS(tmp, 
              file=paste0(rdata_path, sprintf("UCLA_%s.Rds", noise_label)))
      
    } else {
      cat("\nUCLA", noise_proc, "Rds object already exists and --overwrite was not specified. Not writing new Rds object.\n")
    }
  }
}


#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Get breakdown of control vs schizophrenia subjects in cohort

get_dx_breakdown <- function(subject_csv) {
  dx_data <- read.csv(subject_csv) %>%
    mutate(diagnosis = str_to_sentence(diagnosis)) %>%
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

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Plot subject count by diagnosis vs. sex
plot_dx_vs_sex_count <- function(subject_csv) {
  subject_data <- read.csv(subject_csv) %>%
    mutate(diagnosis = str_to_sentence(diagnosis)) %>%
    filter(diagnosis %in% c("Control", "Schz"))
  
  subject_data %>%
    group_by(gender, diagnosis) %>%
    count() %>%
    ggplot(data=., mapping=aes(x=gender, fill=gender, y=n)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=n), vjust=-0.5) +
    facet_grid(diagnosis ~ ., switch="both") +
    scale_y_continuous(expand=c(0,0,0.15,0)) +
    ggtitle("Subject Breakdown by\nSex vs. Diagnosis") +
    ylab("# Subjects") +
    xlab("Sex") +
    theme(legend.position="none",
          plot.title = element_text(hjust=0.5),
          strip.placement = "outside")
}

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Plot age by diagnosis vs. sex
plot_dx_vs_sex_age <- function(subject_csv) {
  subject_data <- read.csv(subject_csv) %>%
    mutate(diagnosis = str_to_sentence(diagnosis)) %>%
    filter(diagnosis %in% c("Control", "Schz"))
  
  subject_data %>%
    ggplot(data=., mapping=aes(x=gender, y=age, fill=gender)) +
    geom_boxplot() +
    facet_grid(diagnosis ~ ., switch="both") +
    ggtitle("Subject Breakdown by\nAge vs. Diagnosis") +
    ylab("Age") +
    xlab("Sex") +
    theme(legend.position="none",
          plot.title = element_text(hjust=0.5),
          strip.placement = "outside")
}

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# catch22 for all regions

catch22_all_regions <- function(TS_df) {
  # Create a new ID that contains both subject ID and brain region
  TS_df$unique_ID <- paste(TS_df$Subject_ID, TS_df$Brain_Region, sep="_")
  
  # Run catch22 using theft
  feature_matrix <- theft::calculate_features(data = TS_df, 
                                              id_var = "unique_ID", 
                                              time_var = "timepoint", 
                                              values_var = "value", 
                                              group_var = "diagnosis", 
                                              feature_set = "catch22",
                                              catch24 = F) %>%
    tidyr::separate(id, into=c("Subject_ID", "Brain_Region"), sep="_")
  
  return(feature_matrix)
}
