Step 3: ROI-Wise catch22 Feature Analysis
================

### Source functions

``` r
source("../helper_functions/visualization_functions.R")
source("ROI_wise_analysis_functions.R")
theme_set(theme_cowplot())
rdata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/Rdata/"
```

### Run multivariable classifier

We will use the multivariable classifier included in `theft` to evaluate
how the catch22 feature set collectively performs at distinguishing
control vs. schizophrenia subjects at each brain region. We will use
z-score normalisation and a linear support vector machine (SVM) with
caret.

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
# If the RDS object doesn't already exist, otherwise load it in
if (!file.exists(paste0(rdata_path, "UCLA_multivar_ROI_res_svmLinear.Rds"))) {
  
  ROI_wise_class_res_list <- list()
  
  for (noise_proc in noise_procs) {
    ROI_wise_class_res <- run_theft_multivar_classifier(rdata_path,
                                                           norm_method = norm_method,
                                                           test_method = test_method,
                                                           noise_proc = noise_proc,
                                                           use_balanced_accuracy = use_balanced_accuracy)
    
    ROI_wise_class_res_list <- rlist::list.append(ROI_wise_class_res_list, ROI_wise_class_res)
  }
  
  ROI_wise_class_res_df <- do.call(plyr::rbind.fill, ROI_wise_class_res_list)
  saveRDS(ROI_wise_class_res_df, file=paste0(rdata_path, "UCLA_multivar_ROI_res_svmLinear.Rds"))
} else {
  ROI_wise_class_res_df <- readRDS(paste0(rdata_path, "UCLA_multivar_ROI_res_svmLinear.Rds"))
}
```

### Examine linear SVM results

We can plot the distribution of accuracy and balanced accuracy values
across the 82 ROIs for each noise processing method:

``` r
# Our three noise-processing methods
noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")
plot_class_acc_w_props(class_res = ROI_wise_class_res_df,
                       rdata_path = rdata_path,
                       noise_procs = noise_procs)
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

This highlights the need for inverse probability weighting so the linear
SVM doesn’t just classify every subject as a control subject to
automatically hit \~70% accuracy.

### In-sample inverse probability weighted SVM with e1071::svm

``` r
noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")

# Use linear SVM
svm_kernel = "linear"

# Run theft's multivariable classifier on each ROI and save to an RDS object
# If the RDS object doesn't already exist, otherwise load it in
if (!file.exists(paste0(rdata_path, "UCLA_e1071_linear_SVM_AROMA_2P.Rds"))) {
  region_wise_SVM_e1071 <- run_e1071_SVM_by_region(rdata_path = rdata_path,
                                                   svm_kernel = svm_kernel,
                                                   noise_procs = noise_procs)
  saveRDS(region_wise_SVM_e1071, file=paste0(rdata_path, "UCLA_e1071_linear_SVM_AROMA_2P.Rds"))
} else {
  region_wise_SVM_e1071 <- readRDS(paste0(rdata_path, "UCLA_e1071_linear_SVM_AROMA_2P.Rds"))
}
```

``` r
# Plot accuracy + balanced accuracy in histograms
# Control subject proportion is highlighted for accuracy, 0.5 is highlighted for balanced accuracy
noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")
plot_class_acc_w_props(class_res = region_wise_SVM_e1071,
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
if (!file.exists(paste0(rdata_path, "UCLA_multivar_ROI_res_svmLinear_inv_prob.Rds"))) {
  
   region_wise_SVM_caret_inv_prob_df <- run_caret_e1071_SVM_by_region(rdata_path = rdata_path,
                                            use_inv_prob_weighting = TRUE,
                                            noise_procs = noise_procs)
  saveRDS(region_wise_SVM_caret_inv_prob_df, file=paste0(rdata_path, "UCLA_multivar_ROI_res_svmLinear_inv_prob.Rds"))
} else {
  region_wise_SVM_caret_inv_prob_df <- readRDS(paste0(rdata_path, "UCLA_multivar_ROI_res_svmLinear_inv_prob.Rds"))
}
```

``` r
plot_class_acc_w_props(class_res = region_wise_SVM_caret_inv_prob_df,
                       rdata_path = rdata_path,
                       noise_procs = noise_procs)
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

As expected with cross-validation, the raw and balanced accuracies are
all shifted lower.