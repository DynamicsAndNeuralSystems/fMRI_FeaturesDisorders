github_dir <- "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/"
rdata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/Rdata/"

# TO-DO: abstract out PCA functions to a helper script
library(tidyverse)
library(sparseSVM)
library(penalizedSVM)
library(cowplot)
theme_set(theme_cowplot())

noise_proc = "AROMA+2P+GMR"
# Clean up names
noise_label <- gsub("\\+", "_", noise_proc)

# Load catch22 data for current noise processing method
feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22.Rds", 
                                                     noise_label)))   

# Prepare data for SVM
data_for_svm <- feature_matrix %>%
  mutate(Unique_ID = paste0(names, "_", Brain_Region),
         .keep = "unused") %>%
  pivot_wider(id_cols = c(Subject_ID, group),
              names_from = Unique_ID, 
              values_from = values) %>%
  drop_na()


svm_mat <- data_for_svm %>%
  dplyr::select(-Subject_ID, -group) %>%
  as.matrix()

################################################################################
# In-sample
################################################################################

# Run L1-regularized in-sample SVM
encoded_group_vec <- ifelse(data_for_svm$group == "Schz", 1, -1)

group_wts <- 1/(prop.table(table(encoded_group_vec)))
in_sample_regularized_svm <- svmfs(x = svm_mat,
                                   y = encoded_group_vec,
                                   fs.method = "1norm",
                                   lambda1.set = c(1),
                                   inner.val.method = "cv",
                                   cross.inner = 10,
                                   grid.search = "discrete",
                                   calc.class.weights = T,
                                   class.weights = group_wts,
                                   parms.coding = "none",
                                   maxIter = 100,
                                   seed = 127)

saveRDS(in_sample_regularized_svm, file=paste0(rdata_path, "Univariate_Combo_Wise_Regularized_SVM_Model.Rds"))

# Evaluate in-sample balanced accuracy
in_sample_preds <- predict(in_sample_regularized_svm, newdata = svm_mat)$pred.class
caret::confusionMatrix(data = in_sample_preds, reference = as.factor(encoded_group_vec))


in_sample_coefs %>%
  group_by(Brain_Region) %>%
  count() %>%
  ungroup() %>%
  mutate(Feature = str_replace_all(Feature, "_", " ")) %>%
  mutate(Feature = fct_reorder(Feature, n)) %>%
  ggplot(data=., mapping = aes(x=Feature, y=n)) +
  geom_bar(aes(fill = Feature), stat="identity") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
  ylab("# Brain Regions") +
  theme(legend.position = "none") +
  coord_flip()


################################################################################
# In-sample
################################################################################

# Specify that group is a factor so that createFolds creates stratified folds
encoded_group_vec <- factor(encoded_group_vec)
group_wts <- 1/(prop.table(table(encoded_group_vec)))

in_sample_regularized_svm <- svmfs(x = svm_mat,
                                   y = encoded_group_vec,
                                   fs.method = "1norm",
                                   lambda1.set = c(1),
                                   inner.val.method = "cv",
                                   cross.inner = 10,
                                   grid.search = "discrete",
                                   calc.class.weights = T,
                                   class.weights = group_wts,
                                   parms.coding = "none",
                                   maxIter = 100,
                                   seed = 127)

saveRDS(in_sample_regularized_svm, file=paste0(rdata_path, "Univariate_Combo_Wise_Regularized_SVM_Model.Rds"))

# Evaluate in-sample balanced accuracy
in_sample_preds <- predict(in_sample_regularized_svm, newdata = svm_mat)$pred.class
caret::confusionMatrix(data = in_sample_preds, reference = as.factor(encoded_group_vec))


in_sample_coefs %>%
  group_by(Brain_Region) %>%
  count() %>%
  ungroup() %>%
  mutate(Feature = str_replace_all(Feature, "_", " ")) %>%
  mutate(Feature = fct_reorder(Feature, n)) %>%
  ggplot(data=., mapping = aes(x=Feature, y=n)) +
  geom_bar(aes(fill = Feature), stat="identity") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
  ylab("# Brain Regions") +
  theme(legend.position = "none") +
  coord_flip()

################################################################################
# Cross-validated
################################################################################

# Create train/test data folds
flds <- caret::createFolds(encoded_group_vec, k = 10, list = TRUE, returnTrain = FALSE)

# Initialize list for performance metrics
in_accuracy_list <- list()
in_balanced_accuracy_list <- list()
out_accuracy_list <- list()
out_balanced_accuracy_list <- list()

k = 10
# Iterate over folds 1 through k
for (i in 1:k) {
  
  # Define test and train data
  test_i <- flds[[i]]
  train_i <- setdiff(1:nrow(svm_mat), test_i)
  
  train_data <- svm_mat[train_i, ]
  train_group <- encoded_group_vec[setdiff(1:nrow(svm_mat), test_i)]
  test_data <- svm_mat[test_i, ]
  test_group <- encoded_group_vec[test_i]
  
  regularized_svm <- svmfs(x = train_data,
                           y = train_group,
                           fs.method = "1norm",
                           lambda1.set = c(1),
                           inner.val.method = "cv",
                           cross.inner = 0,
                           grid.search = "discrete",
                           calc.class.weights = T,
                           class.weights = group_wts,
                           parms.coding = "none",
                           maxIter = 100,
                           seed = 127)
  
  # Generate in-sample predictions based on SVM model
  in_sample_pred <- predict(regularized_svm, newdata = train_data)$pred.class
  train_group <- factor(train_group, levels = levels(in_sample_pred))
  
  # Calculate accuracy and balanced accuracy
  in_accuracy <- sum(in_sample_pred == train_group)/length(in_sample_pred)
  in_balanced_accuracy <- caret::confusionMatrix(reference=train_group, 
                                                 data=in_sample_pred)$byClass[["Balanced Accuracy"]]
  
  # Generate out-of-sample predictions based on SVM model
  out_sample_pred <- predict(regularized_svm, newdata = test_data)$pred.class
  test_group <- factor(test_group, levels = levels(out_sample_pred))
  
  # Calculate accuracy and balanced accuracy
  out_accuracy <- sum(out_sample_pred == test_group)/length(out_sample_pred)
  out_balanced_accuracy <- caret::confusionMatrix(reference=test_group, 
                                                  data=out_sample_pred)$byClass[["Balanced Accuracy"]]
  
  # Append results to list for given fold
  in_accuracy_list <- rlist::list.append(in_accuracy_list, in_accuracy)
  in_balanced_accuracy_list <- rlist::list.append(in_balanced_accuracy_list, in_balanced_accuracy)
  out_accuracy_list <- rlist::list.append(out_accuracy_list, out_accuracy)
  out_balanced_accuracy_list <- rlist::list.append(out_balanced_accuracy_list, out_balanced_accuracy)
  
}



df_res_in_sample <- data.frame(k_fold_iteration = 1:k,
                               accuracy = unlist(in_accuracy_list),
                               balanced_accuracy = unlist(in_balanced_accuracy_list),
                               Sample_Type = "In-sample")
df_res_out_sample <- data.frame(k_fold_iteration = 1:k,
                                accuracy = unlist(out_accuracy_list),
                                balanced_accuracy = unlist(out_balanced_accuracy_list),
                                Sample_Type = "Out-of-sample")

df_res <- plyr::rbind.fill(df_res_in_sample, df_res_out_sample)

df_res %>%
  dplyr::group_by(Sample_Type) %>%
  dplyr::summarise(accuracy_avg = mean(accuracy, na.rm=T),
                   balanced_accuracy_avg = mean(balanced_accuracy, na.rm=T),
                   accuracy_SD = sd(accuracy, na.rm=T),
                   balanced_accuracy_SD = sd(balanced_accuracy, na.rm=T)) %>%
  dplyr::rename("accuracy" = "accuracy_avg",
                "balanced_accuracy" = "balanced_accuracy_avg")