
data_dir <- "/media/sf_Shared_Folder/PhD_work/data/UCLA_Schizophrenia/pydata/AROMA_2P/"
sample_metadata <- "/media/sf_Shared_Folder/PhD_work/data/UCLA_Schizophrenia/participants.csv"
ID_var <- "sampleID"
label_vars <- c("diagnosis")
dim_order <- "ps"
overwrite <- F
yaml_file_base <- "sample.yaml"

npy_files <- list.files(data_dir, pattern="*.npy")
yaml_file <- paste0(data_dir, yaml_file_base)