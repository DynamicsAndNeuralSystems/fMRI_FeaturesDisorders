library(tidyverse)

data_path <- "/headnode1/abry4213/data/HCP100/raw_data/SC_SIFT2_data/"
brain_region_info <- read.csv("/headnode1/abry4213/data/HCP100/Brain_Region_info.csv")

SC_files <- list.files(data_path)

read_SC_data <- function(SC_file, brain_region_info) {
  subject <- gsub("_SC_SIFT2.csv", "", SC_file)
  SC_data <- read.csv(paste0(data_path, SC_file), header=F) 
  colnames(SC_data) <- brain_region_info$Brain_Region
  SC_data$Brain_Region_2 <- brain_region_info$Brain_Region
  SC_data$Sample_ID <- subject
  
  res <- SC_data %>%
    pivot_longer(cols=c(-Sample_ID, -Brain_Region_2),
                 names_to="Brain_Region_1",
                 values_to="nconn") %>%
    filter(Brain_Region_1 != Brain_Region_2)
  
  return(res)
}

full_SC_data <- SC_files %>%
  purrr::map_df(~ read_SC_data(SC_file = .x, brain_region_info = brain_region_info))

# Save SC data to an .Rds file
saveRDS(full_SC_data, file=paste0(data_path, "HCP100_SC_data.Rds"))

# Calculate average # connections between regions
avg_conn_df <- full_SC_data %>%
  group_by(Brain_Region_1, Brain_Region_2) %>%
  summarise(mean_nconn = mean(nconn, na.rm=T)) %>%
  pivot_wider(id_cols="Brain_Region_1", 
              names_from="Brain_Region_2", 
              values_from="mean_nconn") %>%
  column_to_rownames(var="Brain_Region_1")

library(heatmaply)
heatmap(as.matrix(avg_conn_df))
