Step 4: Feature-wise catch22 ROI Analysis
================

### Source functions

``` r
source("../helper_functions/visualization_functions.R")
source("feature_wise_analysis_functions.R")
theme_set(theme_cowplot())
rdata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/Rdata/"
```

### Run multivariable classifier

We will use the multivariable classifier included in `theft` to evaluate
how each catch22 feature performs at distinguishing control
vs. schizophrenia subjects using all 82 brain regions collectively. We
will use z-score normalisation and a linear support vector machine (SVM)
with caret.

The below code chunk calls `run_theft_multivar_classifier`, which is a
wrapper for `fit_multivariable_classifier` from `theft`. Behind the
scenes, we opt to use 10-fold cross validation (`use_k_fold = TRUE`,
`num_folds = 10`) with empirical null model fitting
(`use_empirical_null=TRUE`, `null_testing_method = "null model fits"`)
and gaussian p-value calculation with 10 permutations
(`p_value_method = "gaussian"`, `num_permutations = 10`).

``` r
# Use z-score normalisation
norm_method = "z-score"

# Use linear support vector machine (SVM) classification algorithm
test_method <- "svmLinear"

# Retain balanced accuracy in addition to raw accuracy for each ROI
use_balanced_accuracy <- TRUE

# Run theft's multivariable classifier on each ROI and save to an RDS object
if (!file.exists(paste0(rdata_path, "UCLA_multivar_feature_res_svmLinear.Rds"))) {
  
  feature_wise_class_res_list <- list()
  
  for (noise_proc in noise_procs) {
    feature_wise_class_res <- run_theft_multivar_classifier_by_feature(rdata_path,
                                                                       norm_method = norm_method,
                                                                       test_method = test_method,
                                                                       noise_proc = noise_proc,
                                                                       use_balanced_accuracy = use_balanced_accuracy)
    
    feature_wise_class_res_list <- rlist::list.append(feature_wise_class_res_list, feature_wise_class_res)
  }
  
  feature_wise_class_res_df <- do.call(plyr::rbind.fill, feature_wise_class_res_list)
  saveRDS(feature_wise_class_res_df, file=paste0(rdata_path, "UCLA_multivar_feature_res_svmLinear.Rds"))
} else {
  feature_wise_class_res_df <- readRDS(paste0(rdata_path, "UCLA_multivar_feature_res_svmLinear.Rds"))
}
```

We can plot the distribution of accuracy and balanced accuracy values
across the 22 features for each noise processing method:

``` r
# Our three noise-processing methods
noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")
plot_class_acc_w_props(class_res = feature_wise_class_res_df,
                       rdata_path = rdata_path,
                       noise_procs = noise_procs)
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

### In-sample inverse probability weighted SVM with e1071::svm

``` r
noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")

# Use linear SVM
svm_kernel = "linear"

# Run theft's multivariable classifier on each catch22 feature and save to an RDS object
# If the RDS object doesn't already exist, otherwise load it in
if (!file.exists(paste0(rdata_path, "UCLA_e1071_linear_SVM_by_feature_AROMA_2P.Rds"))) {
  feature_wise_SVM_e1071 <- run_e1071_SVM_by_feature(rdata_path = rdata_path,
                                                   svm_kernel = svm_kernel,
                                                   noise_procs = noise_procs)
  saveRDS(feature_wise_SVM_e1071, file=paste0(rdata_path, "UCLA_e1071_linear_SVM_by_feature_AROMA_2P.Rds"))
} else {
  feature_wise_SVM_e1071 <- readRDS(paste0(rdata_path, "UCLA_e1071_linear_SVM_by_feature_AROMA_2P.Rds"))
}
```

``` r
# Plot accuracy + balanced accuracy in histograms
# Control subject proportion is highlighted for accuracy, 0.5 is highlighted for balanced accuracy
noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")
plot_class_acc_w_props(class_res = feature_wise_SVM_e1071,
                       cv = FALSE,
                       rdata_path = rdata_path,
                       noise_procs = noise_procs)
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

This in-sample inverse probability-weighted SVM with the `e1071` package
performs much better. The distribution of raw accuracies is no longer
largely centered around the control subject proportion point, and the
balanced accuracies are now much higher than 0.5, suggesting that the
classifier is actually attempting to separate schizophrenia vs control
patients.

### 10-fold CV inverse probability weighted SVM with e1071 and caret

``` r
# Try three different noise processing methods
noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")

# Retain balanced accuracy in addition to raw accuracy for each ROI
use_balanced_accuracy <- TRUE

# Implement inverse probability weighting
use_inv_prob_weighting = TRUE

# Run theft's multivariable classifier on each ROI and save to an RDS object
# If the RDS object doesn't already exist, otherwise load it in
if (!file.exists(paste0(rdata_path, "UCLA_multivar_feature_res_svmLinear_inv_prob.Rds"))) {
  
   feature_wise_SVM_caret_inv_prob_df <- run_caret_e1071_SVM_by_feature(rdata_path = rdata_path,
                                            use_inv_prob_weighting = TRUE,
                                            noise_procs = noise_procs)
  saveRDS(feature_wise_SVM_caret_inv_prob_df, file=paste0(rdata_path, "UCLA_multivar_feature_res_svmLinear_inv_prob.Rds"))
} else {
  feature_wise_SVM_caret_inv_prob_df <- readRDS(paste0(rdata_path, "UCLA_multivar_feature_res_svmLinear_inv_prob.Rds"))
}
```

``` r
plot_class_acc_w_props(class_res = feature_wise_SVM_caret_inv_prob_df,
                       rdata_path = rdata_path,
                       noise_procs = noise_procs)
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->