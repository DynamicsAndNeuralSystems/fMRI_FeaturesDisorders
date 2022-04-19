library(theft)
library(tidyverse)


rdata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/Rdata/"

noise_proc = "AROMA+2P"
# Clean up names
noise_label <- gsub("\\+", "_", noise_proc)

# Load catch22 data for current noise processing method
feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22.Rds", 
                                                     noise_label)))      

calculate_balanced_accuracy <- function(data, lev = NULL, model = NULL) {
  # calculate accuracy
  accuracy <- sum(data$pred == data$obs)/length(data$obs)
  
  # Calculate balanced accuracy
  data_cm <- as.data.frame(t(caret::confusionMatrix(data$pred, data$obs)$byClass))
  balanced_accuracy <- mean(data_cm$`Balanced Accuracy`, na.rm=T)
  
  out <- c(accuracy, balanced_accuracy)
  names(out) <- c("Accuracy", "Balanced_Accuracy")
  return(out)
}

# Get control/schz proportions
sample_props <- feature_matrix %>%
  dplyr::group_by(Subject_ID, Brain_Region) %>%
  dplyr::filter(!any(is.na(values))) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(Subject_ID, group) %>%
  dplyr::summarise(control_prop = sum(group=="Control") / n(),
                   schz_prop = sum(group=="Schz")/n())

this_ROI = "ctx-lh-bankssts"
feature_matrix_ROI <- subset(feature_matrix, Brain_Region==this_ROI)

class_res <- fit_multi_feature_classifier(data = feature_matrix_ROI,
                                          id_var = "Subject_ID",
                                          group_var = "group",
                                          by_set = FALSE,
                                          test_method = "svmLinear",
                                          use_balanced_accuracy = TRUE,
                                          use_k_fold = TRUE,
                                          num_folds = 10,
                                          use_empirical_null = TRUE,
                                          null_testing_method = "model free shuffles",
                                          p_value_method = "empirical",
                                          num_permutations = 10000)

# old caret SVM
data_for_svm <- feature_matrix_ROI %>%
  dplyr::rename("id" = "Subject_ID") %>%
  dplyr::select(id, group, names, values) %>%
  tidyr::pivot_wider(id_cols = c("id", "group"), names_from = "names", values_from = "values") %>%
  dplyr::select_if(~sum(!is.na(.)) > 0) %>%
  dplyr::select(where(~dplyr::n_distinct(.) > 1))  %>%
  tidyr::drop_na() %>%
  janitor::clean_names() %>%
  tidyr::pivot_longer(cols = 3:ncol(.), names_to = "names", values_to = "values") %>%
  dplyr::mutate(method = gsub("_.*", "\\1", names)) %>%
  dplyr::mutate(group = make.names(group),
                group = as.factor(group))  %>%
  tidyr::pivot_wider(id_cols = c("id", "group"), names_from = "names", values_from = "values") %>%
  dplyr::select(-c(id))

# Calculate model weights as inverse probability
model_weights <- ifelse(data_for_svm$group == "Control", 
                        1/sample_props$control_prop,
                        1/sample_props$schz_prop)


# Create tune grid
e1071_grid <- expand.grid(cost = 1)
kernlab_grid <- expand.grid(C = 1)

# Train SVM model
fitControl <- caret::trainControl(method = "cv",
                                  number = 10,
                                  savePredictions = "all",
                                  summaryFunction = calculate_balanced_accuracy,
                                  classProbs = TRUE)

# Run e1071 SVM with caret
mod_e1071 <- caret::train(group ~ .,
                    data = data_for_svm,
                    method = "svmLinear2",
                    trControl = fitControl,
                    metric = "Balanced_Accuracy",
                    maximize = T,
                    weights = model_weights,
                    tuneGrid = e1071_grid,
                    preProcess = c("center", "scale", "nzv"))

mod_kernlab <- caret::train(group ~ .,
                            data = data_for_svm,
                            method = "svmLinear",
                            trControl = fitControl,
                            metric = "Balanced_Accuracy",
                            maximize = T,
                            weights = model_weights,
                            tuneGrid = kernlab_grid,
                            preProcess = c("center", "scale", "nzv"))

# Generate predictions
e1071_preds <- mod_e1071$pred
data_for_svm$group <- factor(data_for_svm$group, levels = levels(e1071_preds$pred))
e1071_cm <- list(caret::confusionMatrix(e1071_preds$pred, data_for_svm$group)$table)

kernlab_preds <- mod_kernlab$pred
data_for_svm$group <- factor(data_for_svm$group, levels = levels(kernlab_preds$pred))
kernlab_cm <- list(caret::confusionMatrix(kernlab_preds$pred, data_for_svm$group)$table)

# Get accuracy + balanced accuracy results
e1071_res <- mod_e1071$results %>%
  mutate(Brain_Region = this_ROI,
         Test_Method = "e1071",
         Noise_Proc = noise_proc,
         Confusion_Matrix = I(e1071_cm))

kernlab_res <- mod_kernlab$results %>%
  mutate(Brain_Region = this_ROI,
         Test_Method = "kernlab",
         Noise_Proc = noise_proc,
         Confusion_Matrix = I(kernlab_cm))


plyr::rbind.fill(e1071_res, kernlab_res)
