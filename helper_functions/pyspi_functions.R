read_pyspi_pkl_into_RDS <- function(data_path,
                                    noise_procs = c("AROMA+2P",
                                                    "AROMA+2P+GMR",
                                                    "AROMA+2P+DiCER")) {
  # Iterate over each noise-processing method
  for (noise_proc in noise_procs) {
    noise_label = gsub("\\+", "_", noise_proc)
    
    # Define data path for this noise processing method
    np_data_path <- paste0(data_path, noise_label, "/")
    
    # Initialize list to store data
    subj_data_list <- list()
    
    # Iterate over each subject
    for (subject in unique(list.dirs(np_data_path, recursive = F, full.names = F))) {
      
      # If subject doesn't have a corresponding pyspi RDS file for this 
      # noise-processing method, create one
      if (!file.exists(paste0(np_data_path, subject, "_pyspi.Rds"))) {
        tryCatch({subject_pkl_data <- extract_df_from_pkl(paste0(np_data_path, subject, "/calc.pkl")) %>%
          mutate(Subject_ID = subject,
                 group = subset(subject_csv, sampleID == subject) %>% pull(diagnosis))
        
        # Save results to an RDS file for this subject
        saveRDS(subject_pkl_data, file=paste0(np_data_path, subject, "_pyspi.Rds"))},
        error = function(e){
          cat("\nError for subject", subject, "\n")
          print(e)
        })
        # 
        # # Append results to list
        # subj_data_list <- rlist::list.append(subj_data_list, subject_pkl_data)
      }
    }
    
    # # Merge subjects' data together
    # all_pyspi_data <- do.call(plyr::rbind.fill, subj_data_list)
    # saveRDS(all_pyspi_data, paste0(data_path, "All_subject_pyspi_", noise_label, ".Rds"))
  }
}