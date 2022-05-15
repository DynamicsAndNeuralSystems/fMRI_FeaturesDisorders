# Load packages
library(purrr)
library(dplyr)
library(stringr)

# Data paths
data_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/"
rdata_path <- paste0(data_path, "Rdata/")

# Use the three noise processing methods
noise_procs <- c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")

# Define output path for CSV files
CSV_input_path <- paste0(data_path, "catchaMouse16_data/output_CSV/")


# Iterate over each noise-processing method
read_in_noise_proc_data <- function(noise_proc) {
  noise_label = gsub("\\+", "_", noise_proc)
  if (!file.exists(paste0(rdata_path, "UCLA_",
                          noise_label, "_catchaMouse16.Rds"))) {
    np_files <- grep(paste0(noise_label, "_cm16.csv"), 
                     list.files(CSV_input_path), value = TRUE)
    
    noise_proc_list <- list()
    for (file in np_files) {
      subject <- str_split(file, "_", n=2)[[1]][1]
      brain_region <- str_split(file, "_")[[1]][2]
      noise_proc <- noise_proc
      
      df <- read.csv(paste0(CSV_input_path, file),
                     header=F, na.strings = "NaN")[,1:2]
      colnames(df) <- c("Value",
                        "catchaMouse16_Feature")
      
      df$catchaMouse16_Feature <- gsub("[[:space:]]", "", df$catchaMouse16_Feature)
      
      df <- df %>%
        mutate(Subject_ID = subject,
               Brain_Region = brain_region,
               Noise_Proc = noise_proc)
      
      noise_proc_list <- rlist::list.append(noise_proc_list, df)
    }
    noise_proc_res <- do.call(plyr::rbind.fill,
                              noise_proc_list)
    
    saveRDS(noise_proc_res, file=paste0(rdata_path, "UCLA_",
                                        noise_label, "_catchaMouse16.Rds"))
  }
}

noise_procs %>%
  purrr::map_df( ~ read_in_noise_proc_data(noise_proc = .x))



