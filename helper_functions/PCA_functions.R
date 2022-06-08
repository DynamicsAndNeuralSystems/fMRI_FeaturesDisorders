#------------------------------------
# This script compiles movement data across the UCLA subjects
#------------------------------------

#--------------------------------------
# Author: Annie Bryant, 5 June 2022
#--------------------------------------

library(FactoMineR)
library(factoextra)

#-------------------------------------------------------------------------------
# Function to read in univariate TS feature data and return subjects with NA 
# values. For each noise-processing method, the number of TS features with all 
# NA for all brain regions are given.
#-------------------------------------------------------------------------------

run_univ_PCA_by_group_var <- function(feature_matrix,
                                      grouping_variable = "Brain_Region",
                                      feature_var = "names") {
  PCA_object_list <- list()
  PCA_scores_list <- list()
  PCA_eigenvalues_list <- list()
  
  grouping_var_vector <- feature_matrix %>%
    dplyr::pull(get(grouping_variable)) %>%
    unique()
  
  for (this_group in grouping_var_vector) {
    df_for_PCA <- subset(feature_matrix, 
                         get(grouping_variable) == this_group) %>%
      pivot_wider(id_cols = c(Subject_ID, group, 
                              grouping_variable),
                  names_from = get(feature_var),
                  values_from = values) %>%
      drop_na()
  }
}


for (this_ROI in unique(feature_matrix$Brain_Region)) {
  ROI_data <- subset(feature_matrix, Brain_Region == this_ROI) %>%
    pivot_wider(id_cols = c(Subject_ID, group, Brain_Region),
                names_from = names,
                values_from = values) %>%
    drop_na()
  ROI_data_mat <- ROI_data %>%
    select(-Subject_ID, -group, -Brain_Region) %>%
    as.matrix()
  
  ROI_data_PCA <- prcomp(ROI_data_mat, center = TRUE, scale. = TRUE)
  ROI_PCA_list[[this_ROI]] <- ROI_data_PCA
  
  ROI_PC_vals <- as.data.frame(ROI_data_PCA$x)
  ROI_PC_vals$group <- ROI_data$group
  ROI_PC_vals$Brain_Region <- this_ROI
  ROI_PC_scores_list <- rlist::list.append(ROI_PC_scores_list, ROI_PC_vals)
  
  ROI_eigen  <- get_eig(ROI_data_PCA) %>%
    mutate(Components = 1:nrow(.)) %>%
    dplyr::rename("Percent_Variance" = "variance.percent",
                  "Cumulative_Variance" = "cumulative.variance.percent") %>%
    mutate(Brain_Region = this_ROI)
  
  ROI_eigen_list <- rlist::list.append(ROI_eigen_list, ROI_eigen)
}

ROI_eigen_df <- do.call(plyr::rbind.fill, ROI_eigen_list)
ROI_PC_scores <- do.call(plyr::rbind.fill, ROI_PC_scores_list)