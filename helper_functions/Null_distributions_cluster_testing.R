# Parse arguments

pairwise_data_file="/project/hctsa/annie/data/scz/UCLA/Rdata/pyspi_SPI_pairwise_CV_linear_SVM_pyspi_19_inv_prob.Rds"
SPI_directionality_file="/project/hctsa/annie/github/fMRI_FeaturesDisorders/pairwise_analysis/SPI_Direction_Info.csv"
rdata_path="/project/hctsa/annie/data/scz/UCLA/Rdata/"
github_dir="/project/hctsa/annie/github/fMRI_FeaturesDisorders/"
null_iter_number=1
feature_set="pyspi_19"
svm_kernel="linear"
grouping_var="SPI"
svm_feature_var="region_pair"
test_package="e1071"
noise_proc="AROMA+2P+GMR"
return_all_fold_metrics=TRUE
use_inv_prob_weighting=TRUE
use_SMOTE=FALSE

# Load data
cat("\nNow reading in pyspi RDS file\n")
pairwise_data <- readRDS(pairwise_data_file)
cat("\nNow reading in SPI directionality file\n")
SPI_directionality <- read.csv(SPI_directionality_file)

# Source linear SVM functions
source(paste0(github_dir, "helper_functions/Linear_SVM.R"))

# Define output directory
if (use_inv_prob_weighting) {
  output_dir <- paste0(rdata_path, sprintf("Pairwise_%s_inv_prob_null_model_fits/",
                                           feature_set))
} else {
  output_dir <- paste0(rdata_path, sprintf("Pairwise_%s_unweighted_null_model_fits/",
                                           feature_set))
}

icesTAF::mkdir(output_dir)


# Run null iteration
cat("\nNow running null iteration:\n")
# null_out <- run_pairwise_cv_svm_by_input_var(pairwise_data = pairwise_data,
#                                              SPI_directionality = SPI_directionality,
#                                              svm_kernel = svm_kernel,
#                                              grouping_var = grouping_var,
#                                              svm_feature_var = svm_feature_var,
#                                              test_package = test_package,
#                                              noise_proc = noise_proc,
#                                              return_all_fold_metrics = return_all_fold_metrics,
#                                              use_inv_prob_weighting = use_inv_prob_weighting,
#                                              use_SMOTE = use_SMOTE,
#                                              shuffle_labels = TRUE)
svm_feature_var_name = svm_feature_var
grouping_var_name = "SPI"
grouping_var_vector <- unique(pairwise_data$SPI)

cat("\nHead of pairwise data:\n")
head(pairwise_data)

cat("\n\nHead of SPI directionality data:\n")
head(SPI_directionality)

# Filter by directionality
pairwise_data <- pairwise_data %>%
  left_join(., SPI_directionality) %>%
  rowwise() %>%
  mutate(region_pair = case_when(Direction == "Undirected" ~ ifelse(brain_region_1 < brain_region_2,
                                                                    paste0(brain_region_1, "_", brain_region_2),
                                                                    paste0(brain_region_2, "_", brain_region_1)),
                                 Direction == "Directed" ~ paste0(brain_region_1, "_", brain_region_2))) %>%
  dplyr::select(-brain_region_1, -brain_region_2)  %>%
  distinct(Subject_ID, SPI, region_pair, .keep_all = T)


for (group_var in unique(grouping_var_vector)) {
  if (grouping_var == "Combo") {
    data_for_SVM <- pairwise_data %>%
      dplyr::select(Subject_ID, group, Combo, value) %>%
      tidyr::pivot_wider(id_cols = c(Subject_ID, group),
                         names_from = Combo,
                         values_from 
                         = value) %>%
      dplyr::select(-Subject_ID) %>%
      # Drop columns that are all NA/NAN
      dplyr::select(where(function(x) any(!is.na(x)))) %>%
      # Drop rows with NA for one or more column
      drop_na()
    
  } else {
    # Otherwise iterate over each separate group
    data_for_SVM <- subset(pairwise_data, get(grouping_var_name) == group_var) %>%
      dplyr::ungroup() %>%
      dplyr::select(Subject_ID, group, svm_feature_var_name, value) %>%
      tidyr::pivot_wider(id_cols = c(Subject_ID, group),
                         names_from = svm_feature_var_name,
                         values_from 
                         = value) %>%
      dplyr::select(-Subject_ID) %>%
      # Drop columns that are all NA/NAN
      dplyr::select(where(function(x) any(!is.na(x)))) %>%
      # Drop rows with NA for one or more column
      drop_na()
  }
  
  # Define sample weights
  # Default is 1 and 1 if use_inv_prob_weighting is not included
  if (use_inv_prob_weighting) {
    # Get control/schz proportions
    sample_wts <- as.list(1/prop.table(table(data_for_SVM$group)))
  } else {
    sample_wts <- list("Control" = 1, "Schz" = 1)
  }
  
  if (nrow(data_for_SVM) > 0) {
    # Run k-fold linear SVM
    SVM_results <- k_fold_CV_linear_SVM(input_data = data_for_SVM,
                                        k = 10,
                                        svm_kernel = svm_kernel,
                                        sample_wts = sample_wts,
                                        use_SMOTE = use_SMOTE,
                                        shuffle_labels = TRUE,
                                        return_all_fold_metrics = return_all_fold_metrics) %>%
      dplyr::mutate(grouping_var = group_var,
                    Noise_Proc = noise_proc,
                    use_inv_prob_weighting = use_inv_prob_weighting,
                    use_SMOTE = use_SMOTE)
    
    # Append results to list
    class_res_list <- rlist::list.append(class_res_list,
                                         SVM_results)
  } else {
    cat("\nNo observations available for", group_var, "after filtering.\n")
  }
}

# Combine results from all regions into a dataframe
null_out <- do.call(plyr::rbind.fill, class_res_list)

# Save null results to RDS
if (use_inv_prob_weighting) {
  saveRDS(null_out, file=sprintf("%s/Pairwise_%s_inv_prob_null_model_fit_iter_%s.Rds",
                                 output_dir, feature_set, null_iter_number))
} else {
  saveRDS(null_out, file=sprintf("%s/Pairwise_%s_unweighted_null_model_fit_iter_%s.Rds",
                                 output_dir, feature_set, null_iter_number))
}
