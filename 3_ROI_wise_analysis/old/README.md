Region-Wise Classification Analysis
================

### Source functions

``` r
source("../helper_functions/visualization_functions.R")
source("region_by_region_analysis.R")
rdata_path <- "D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/Rdata/"
```

We can visualize PC1 vs. PC2 for an example ROI across the three
different noise processing methods, like the left isthmus cingulate:

``` r
# Set an example ROI
this_ROI = "ctx-lh-isthmuscingulate"
ROI_label = "Left_isthmus_cingulate"

# Plot PC1 vs PC2 for the left entorhinal cortex as an example
PCA_dimplot_for_ROI(this_ROI = this_ROI,
                    rdata_path = rdata_path)
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

We can also visualize the violin plots for the raw data and the two
normalisation methods:

``` r
# Set an example ROI
this_ROI = "ctx-lh-isthmuscingulate"

# Just use AROMA+2P for example
noise_proc = "AROMA+2P"
noise_label <- gsub("\\+", "_", noise_proc)
noise_label <- gsub("\\+", "_", noise_proc)
violin_plots_for_ROI(this_ROI = this_ROI,
                     rdata_path = rdata_path,
                     noise_proc = noise_proc)
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

### Manual t-test

Let’s start with a very simple t-test for catch22 feature values in
control vs schizophrenia subjects by brain region:

``` r
t_test_res <- t_test_by_region(rdata_path=rdata_path) %>%
  dplyr::rename("statistic_value"="estimate",
                "feature"="names") 

# Just use AROMA+2P for example
noise_proc = "AROMA+2P"
noise_label <- gsub("\\+", "_", noise_proc)

t_stat_histograms(t_test_res=t_test_res,
                  noise_proc = noise_proc)
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Let’s take a birds-eye view of the T-test classification metrics in with
each of the noise-processing datasets.

``` r
statistic = "T statistic"
test_method = "t_test"

norm_methods <- c("non-normalised", "z-score", "RobustSigmoid")

min_t <- min(t_test_res$statistic_value, na.rm=T)
max_t <- max(t_test_res$statistic_value, na.rm=T)

if (abs(min_t) > max_t) {
  min_val = min_t
  max_val = abs(min_t)
} else {
  min_val = -1*max_t
  max_val = max_t
}

# Just use AROMA+2P for this example
noise_proc <- "AROMA+2P"
noise_label <- gsub("\\+", "_", noise_proc)

for (norm_method in norm_methods) {
  classification_results <- t_test_res %>%
    filter(Noise_Proc == noise_proc,
           Norm_Method == norm_method)
  
  # Heatmap
  plot_class_stat_heatmap(classification_results = classification_results,
                          statistic = statistic,
                          noise_proc = noise_proc,
                          norm_method = norm_method,
                          min_val = min_val,
                          max_val = max_val)
}
```

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->

### Theft

Run univariate classification on each brain region separately, using
both t-test and linear SVM:

``` r
# Run univariate classification on each brain region separately
test_methods <- c("t-test", "svmLinear")
norm_methods <- c("z-score", "RobustSigmoid")


for (norm_method in norm_methods) {
  for (test_method in test_methods) {
    for (noise_proc in c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER")) { 
      
      # Clean up names
      test_label <- gsub("-", "_", test_method)
      noise_label <- gsub("\\+", "_", noise_proc)
      
      # Only run if RDS file doesn't yet exist
      if (!file.exists(paste0(rdata_path, "UCLA_", 
                              noise_label, "_catch22_ROIwise_",  
                              norm_method, "_",
                              test_label, ".Rds"))) {
        
        # Load catch22 feature matrix
        feature_matrix <- readRDS(paste0(rdata_path, sprintf("UCLA_%s_catch22.Rds", 
                                                             noise_label)))
        
        # Run ROI-by-ROI analysis
        run_region_by_region_analysis(feature_matrix = feature_matrix, 
                                      test_method = test_method, 
                                      norm_method = norm_method,
                                      rdata_path = rdata_path, 
                                      noise_proc = noise_proc)
      }
      # clean up memory
      gc()
    }
    
  }
}
```

Let’s take a birds-eye view of the T-test classification metrics in with
each of the noise-processing datasets.

``` r
test_statistic = "T statistic"
test_method = "t_test"

norm_methods <- c("z-score", "RobustSigmoid")

class_res_list <- list()
for (noise_proc in c("AROMA+2P",
                     "AROMA+2P+GMR",
                     "AROMA+2P+DiCER")) {
  noise_label <- gsub("\\+", "_", noise_proc)
  
  for (norm_method in norm_methods) {
    classification_results <- readRDS(paste0(rdata_path, 
                                             sprintf("UCLA_%s_catch22_ROIwise_%s_%s.Rds",
                                                     noise_label,
                                                     norm_method,
                                                     test_method)))
    
    # Heatmap
    png(sprintf("plots/UCLA_%s_catch22_ROIwise_%s_%s.png",
                noise_label, test_method, norm), width=12, height=7, units="in", res=300)
    plot_class_stat_heatmap(classification_results = classification_results,
                            statistic = test_statistic,
                            noise_proc = noise_proc,
                            norm_method = norm)
    dev.off()
  }
}
```
