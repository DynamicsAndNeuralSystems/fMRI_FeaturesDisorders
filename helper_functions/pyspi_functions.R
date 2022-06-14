#-------------------------------------------------------------------------------
# Function to read in pySPI calc.pkl files per subject and convert
# The SPI result data into an RDS file per subject
#-------------------------------------------------------------------------------
read_pyspi_pkl_into_RDS <- function(data_path,
                                    subject_csv_file,
                                    noise_procs = c("AROMA+2P",
                                                    "AROMA+2P+GMR",
                                                    "AROMA+2P+DiCER")) {
  
  # Load subject data
  subject_data <- read.csv(subject_csv_file)
  # Iterate over each noise-processing method
  for (noise_proc in noise_procs) {
    noise_label = gsub("\\+", "_", noise_proc)
    
    # Define data path for this noise processing method
    np_data_path <- paste0(data_path, noise_label, "/")
    
    # Iterate over each subject
    for (subject in unique(list.dirs(np_data_path, recursive = F, full.names = F))) {
      
      # If subject doesn't have a corresponding pyspi RDS file for this 
      # noise-processing method, create one
      if (!file.exists(paste0(np_data_path, subject, "_pyspi.Rds"))) {
        cat("\nNow prepping data for", subject, noise_proc, "\n")
        tryCatch({subject_pkl_data <- extract_df_from_pkl(paste0(np_data_path, subject, "/calc.pkl")) %>%
          mutate(Subject_ID = subject,
                 group = subset(subject_data, sampleID == subject) %>% pull(diagnosis))

        # Save results to an RDS file for this subject
        saveRDS(subject_pkl_data, file=paste0(np_data_path, subject, "_pyspi.Rds"))},
        error = function(e){
          cat("\nError for subject", subject, "\n")
          print(e)
        })
      }
    }
  }
}

#-------------------------------------------------------------------------------
# Function to read in pySPI calc.pkl files per subject and convert
# The SPI result data into an RDS file per subject
#-------------------------------------------------------------------------------
merge_pyspi_res_for_study <- function(data_path,
                                      input_dataset_name = "UCLA",
                                      ROI_index_file,
                                      noise_procs = c("AROMA+2P",
                                                      "AROMA+2P+GMR",
                                                      "AROMA+2P+DiCER")) {
  # Read in ROI index data
  ROI_index <- read.csv(ROI_index_file)
  
  # Iterate over each noise-processing method
  for (noise_proc in noise_procs) {
    noise_label = gsub("\\+", "_", noise_proc)
    if (!file.exists(paste0(data_path, input_dataset_name, 
                            "_all_subject_pyspi_", noise_label, ".Rds"))) {
      # Define data path for this noise processing method
      np_data_path <- paste0(data_path, noise_label, "/")
      
      # Initialize list to store data
      subj_data_list <- list()
      
      # Iterate over each subject
      subjects = unique(list.dirs(np_data_path, recursive = F, full.names = F))
      for (subject in subjects) {
        
        # If subject doesn't have a corresponding pyspi RDS file for this 
        # noise-processing method, create one
        tryCatch({
          subject_pyspi_res <- readRDS(paste0(np_data_path, subject, "_pyspi.Rds")) 
          
          # Append results to list
          subj_data_list <- rlist::list.append(subj_data_list, subject_pyspi_res)
        },
        error = function(e){
          cat("\nError for subject", subject, "\n")
          print(e)
        })
        
      }
      
      # Merge subjects' data together
      all_pyspi_data <- do.call(plyr::rbind.fill, subj_data_list)  %>%
        mutate(comparison = row_number(),
               group = stringr::str_to_sentence(group)) %>%
        pivot_longer(cols = c(brain_region_1,
                              brain_region_2),
                     names_to = "Region_Number",
                     values_to = "Index") %>%
        left_join(ROI_index) %>%
        dplyr::select(-Index) %>%
        pivot_wider(id_cols = c("Subject_ID", "group", "SPI", "value", "comparison"),
                    names_from = "Region_Number",
                    values_from = "ROI") %>%
        dplyr::select(-comparison) 
      
      saveRDS(all_pyspi_data, paste0(data_path, input_dataset_name, 
                                     "_all_subject_pyspi_", noise_label, ".Rds"))
    }
  }
}