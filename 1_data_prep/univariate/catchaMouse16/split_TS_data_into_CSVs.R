# Load dplyr
library(purrr)
library(dplyr)
# Data paths
data_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/"
rdata_path <- paste0(data_path, "Rdata/")

# Use the three noise processing methods
noise_procs <- c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")

# Define output path for CSV files
CSV_output_path <- paste0(data_path, "catchaMouse16_data/input_CSV/")

# Define function to write TS data to a CSV per subject
write_subject_ROI_csv <- function(unique_ID, TS_data, noise_label) {
  subject_ROI_data <- subset(TS_data, Unique_ID==unique_ID) %>%
    pull(value)
  write.table(subject_ROI_data, 
              file=paste0(CSV_output_path, unique_ID, "_", noise_label, "_TS.csv"),
              sep=",",
              row.names = FALSE,
              col.names = FALSE)
}

for (noise_proc in noise_procs) {
  noise_label = gsub("\\+", "_", noise_proc)
  
  TS_df <- readRDS(paste0(rdata_path, "UCLA_", noise_label, ".Rds"))
  
  TS_df <- TS_df %>%
    dplyr::select(Subject_ID, Brain_Region, value) %>%
    tidyr::unite("Unique_ID", c(Subject_ID, Brain_Region))
  
  unique(TS_df$Unique_ID) %>%
    purrr::map(~ write_subject_ROI_csv(unique_ID = .x,
                                   TS_data = TS_df,
                                   noise_label = noise_label))
  
}