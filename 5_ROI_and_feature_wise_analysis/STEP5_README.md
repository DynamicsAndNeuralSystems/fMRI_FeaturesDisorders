Step 5: ROI+Feature Classification Analysis
================

### Source functions

## Cross-validated SVM classification

### 10-fold cross-validated linear SVM

We can implement 10-fold cross-validation (CV) with the `caret` package.

![](Step5_README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

As with in-sample SVM, the unweighted input samples are virtually all
classified as control subjects across all 82 ROIs using the 10-fold
cross-validation linear SVM with caret.

### 10-fold cross-validated linear SVM with inverse probability weighting

![](Step5_README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Surprisingly, incorporating inverse probability weighting has minimal
impact when it comes to the ten-fold cross-validated SVM. Of note, the
in-sample and cross-validated SVM were both run with kernlab::ksvm using
default parameters.

### 10-fold cross-validated linear SVM with SMOTE

![](Step5_README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

Surprisingly, incorporating inverse probability weighting has minimal
impact when it comes to the ten-fold cross-validated SVM. Of note, the
in-sample and cross-validated SVM were both run with kernlab::ksvm using
default parameters.

## Model-free shuffle null distribution

### Generating null distributions from model-free shuffles

This first model-free shuffles method is borrowed from Trent’s
implementation in theft. With this method, the input class labels (Schz
or Control) are randomly shuffled N times, and for each iteration, the
classification accuracy and balanced accuracy are calculated. This
yields a null distribution of accuracies and balanced accuracies,
circumventing the need for running any classification algorithms across
iterations.

Here, I’ve run 1,000,000 iterations of the model-free shuffle,
generating 1,000,000 null values for Accuracy and Balanced Accuracy,
respectively. Since this method is independent of ROI/feature combo, the
same null distribution can be used to compare with each ROI/feature
combo separately.

### CV linear SVM

I’ve plotted the distribution of null accuracies (gray) alongside the
actual accuracies (green, red, and blue) for the ROI/Feature combos per
noise-processing method.

![](Step5_README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Noise_Proc
</th>
<th style="text-align:right;">
balanced_accuracy
</th>
<th style="text-align:right;">
bal_acc_p
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
AROMA+2P
</td>
<td style="text-align:right;font-weight: bold;color: black !important;background-color: palegreen !important;">
0.6256061
</td>
<td style="text-align:right;font-weight: bold;color: black !important;background-color: palegreen !important;">
0.000676
</td>
</tr>
<tr>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
AROMA+2P+GMR
</td>
<td style="text-align:right;font-weight: bold;color: black !important;background-color: palegreen !important;">
0.6339394
</td>
<td style="text-align:right;font-weight: bold;color: black !important;background-color: palegreen !important;">
0.000167
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+DiCER
</td>
<td style="text-align:right;">
0.5093182
</td>
<td style="text-align:right;">
0.368329
</td>
</tr>
</tbody>
</table>

This table summarises the number of ROIs for which raw accuracy or
balanced accuracy is significantly greater than the model-free shuffle
null distribution, both before and after adjusting for multiple
comparisons with BH-FDR.

### CV linear SVM – inv prob

![](Step5_README_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

I’ve plotted the distribution of null accuracies (gray) alongside the
actual accuracies (green, red, and blue) for the ROI/Feature combos per
noise-processing method.

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Noise_Proc
</th>
<th style="text-align:right;">
balanced_accuracy
</th>
<th style="text-align:right;">
bal_acc_p
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
AROMA+2P
</td>
<td style="text-align:right;font-weight: bold;color: black !important;background-color: palegreen !important;">
0.6088636
</td>
<td style="text-align:right;font-weight: bold;color: black !important;background-color: palegreen !important;">
0.002201
</td>
</tr>
<tr>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
AROMA+2P+GMR
</td>
<td style="text-align:right;font-weight: bold;color: black !important;background-color: palegreen !important;">
0.6389394
</td>
<td style="text-align:right;font-weight: bold;color: black !important;background-color: palegreen !important;">
0.000167
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+DiCER
</td>
<td style="text-align:right;">
0.5434848
</td>
<td style="text-align:right;">
0.143590
</td>
</tr>
</tbody>
</table>

This table summarises the number of ROIs for which raw accuracy or
balanced accuracy is significantly greater than the model-free shuffle
null distribution, both before and after adjusting for multiple
comparisons with BH-FDR.

### CV linear SVM – SMOTE

![](Step5_README_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

I’ve plotted the distribution of null accuracies (gray) alongside the
actual accuracies (green, red, and blue) for the ROI/Feature combos per
noise-processing method.

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Noise_Proc
</th>
<th style="text-align:right;">
balanced_accuracy
</th>
<th style="text-align:right;">
bal_acc_p
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
AROMA+2P
</td>
<td style="text-align:right;font-weight: bold;color: black !important;background-color: palegreen !important;">
0.5772727
</td>
<td style="text-align:right;font-weight: bold;color: black !important;background-color: palegreen !important;">
0.016072
</td>
</tr>
<tr>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
AROMA+2P+GMR
</td>
<td style="text-align:right;font-weight: bold;color: black !important;background-color: palegreen !important;">
0.6296212
</td>
<td style="text-align:right;font-weight: bold;color: black !important;background-color: palegreen !important;">
0.000636
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+DiCER
</td>
<td style="text-align:right;">
0.5205303
</td>
<td style="text-align:right;">
0.241621
</td>
</tr>
</tbody>
</table>

This table summarises the number of ROIs for which raw accuracy or
balanced accuracy is significantly greater than the model-free shuffle
null distribution, both before and after adjusting for multiple
comparisons with BH-FDR.

## Empirical model-based pooled null distribution

### Generating null distributions from pooled null model fits

In contrast to the model-free shuffle method, here we are actually
shuffling the input class labels right before running the linear SVM
over N=100 iterations per ROI (N=82) and pooling the resulting accuracy
and balanced accuracy values, to generate empirical null distributions
of N=1,000 data points each, respectively.

### Unweighted

### CV, inverse probability weighted

![](Step5_README_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

The fitted empirical null model distribution is fairly similar to the
real accuracy and balanced accuracy values using in-sample linear SVM
with no reweighting.

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Noise_Proc
</th>
<th style="text-align:right;">
balanced_accuracy
</th>
<th style="text-align:left;">
bal_acc_p
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
AROMA+2P
</td>
<td style="text-align:right;font-weight: bold;color: black !important;background-color: palegreen !important;">
1.0000000
</td>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
0e+00
</td>
</tr>
<tr>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
AROMA+2P
</td>
<td style="text-align:right;font-weight: bold;color: black !important;background-color: palegreen !important;">
0.6088636
</td>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
1e-02
</td>
</tr>
<tr>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
AROMA+2P+GMR
</td>
<td style="text-align:right;font-weight: bold;color: black !important;background-color: palegreen !important;">
1.0000000
</td>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
0e+00
</td>
</tr>
<tr>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
AROMA+2P+GMR
</td>
<td style="text-align:right;font-weight: bold;color: black !important;background-color: palegreen !important;">
0.6389394
</td>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
0e+00
</td>
</tr>
<tr>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
AROMA+2P+DiCER
</td>
<td style="text-align:right;font-weight: bold;color: black !important;background-color: palegreen !important;">
1.0000000
</td>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
0.0e+00
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+DiCER
</td>
<td style="text-align:right;">
0.5434848
</td>
<td style="text-align:left;">
1.3e-01
</td>
</tr>
</tbody>
</table>

### CV, SMOTE

![](Step5_README_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

The fitted empirical null model distribution is fairly similar to the
real accuracy and balanced accuracy values using in-sample linear SVM
with no reweighting.

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Noise_Proc
</th>
<th style="text-align:right;">
balanced_accuracy
</th>
<th style="text-align:left;">
bal_acc_p
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
AROMA+2P
</td>
<td style="text-align:right;font-weight: bold;color: black !important;background-color: palegreen !important;">
1.0000000
</td>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
0e+00
</td>
</tr>
<tr>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
AROMA+2P
</td>
<td style="text-align:right;font-weight: bold;color: black !important;background-color: palegreen !important;">
0.5772727
</td>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
1e-02
</td>
</tr>
<tr>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
AROMA+2P+GMR
</td>
<td style="text-align:right;font-weight: bold;color: black !important;background-color: palegreen !important;">
1.0000000
</td>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
0e+00
</td>
</tr>
<tr>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
AROMA+2P+GMR
</td>
<td style="text-align:right;font-weight: bold;color: black !important;background-color: palegreen !important;">
0.6296212
</td>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
0e+00
</td>
</tr>
<tr>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
AROMA+2P+DiCER
</td>
<td style="text-align:right;font-weight: bold;color: black !important;background-color: palegreen !important;">
1.0000000
</td>
<td style="text-align:left;font-weight: bold;color: black !important;background-color: palegreen !important;">
0.0e+00
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+DiCER
</td>
<td style="text-align:right;">
0.5205303
</td>
<td style="text-align:left;">
1.9e-01
</td>
</tr>
</tbody>
</table>

## Comparing model-free shuffle with pooled empirical null distributions

![](Step5_README_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->
