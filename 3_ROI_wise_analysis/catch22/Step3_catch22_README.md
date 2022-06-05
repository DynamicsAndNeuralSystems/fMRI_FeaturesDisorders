Step 3: ROI-Wise catch22 Feature Analysis
================

### Source functions

## Cross-validated SVM classification

### 10-fold cross-validated linear SVM

I have chosen to use 10-fold cross validation via manual implementation,
as the sample reweighting options in caret were limited and difficult to
interpret.

![](Step3_catch22_README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

Interestingly, unlike the in-sample results, there is a fair spread of
accuracy and balanced accuracy values outside of the proportions
expected from classifying all subjects as controls.

However, there still is a balanced accuracy peak around 0.5, so we move
forward with inverse probability reweighting.

### 10-fold cross-validated linear SVM with inverse probability weighting

![](Step3_catch22_README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

As with the in-sample results, the accuracy values are negatively
shifted while the balanced accuracy values are positively shifted after
applying inverse probability reweighting to the samples.

### 10-fold cross-validated linear SVM with SMOTE

sa\[’\]

![](Step3_catch22_README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

As with the in-sample results, the accuracy values are negatively
shifted while the balanced accuracy values are positively shifted after
applying inverse probability reweighting to the samples.

## Model-free shuffle null distribution

### Generating null distributions from model-free shuffles

This first model-free shuffles method is borrowed from Trent’s
implementation in theft. With this method, the input class labels (Schz
or Control) are randomly shuffled N times, and for each iteration, the
classification accuracy and balanced accuracy are calculated. This
yields a null distribution of accuracies and balanced accuracies,
circumventing the need for running any classification algorithms across
iterations.

Here, I’ve run 100,000 iterations of the model-free shuffle, generating
100,000 null values for Accuracy and Balanced Accuracy, respectively.
Since this method is independent of brain region, the same null
distribution can be used to compare with each brain region separately.

### Unweighted linear SVM

![](Step3_catch22_README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

I’ve plotted the distribution of null accuracies (teal) alongside the
actual accuracies (pink) for the 82 ROIs on the left. Let’s zoom in on
AROMA+2P and pick the five brain regions with the highest
cross-validated balanced accuracy:

![](Step3_catch22_README_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Noise_Proc
</th>
<th style="text-align:left;">
Sample_Type
</th>
<th style="text-align:right;">
num_sig_acc
</th>
<th style="text-align:right;">
num_sig_acc_fdr
</th>
<th style="text-align:right;">
num_sig_bacc
</th>
<th style="text-align:right;">
num_sig_bacc_fdr
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
AROMA+2P
</td>
<td style="text-align:left;">
In-sample
</td>
<td style="text-align:right;">
82
</td>
<td style="text-align:right;">
82
</td>
<td style="text-align:right;">
51
</td>
<td style="text-align:right;">
45
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P
</td>
<td style="text-align:left;">
Out-of-sample
</td>
<td style="text-align:right;">
68
</td>
<td style="text-align:right;">
68
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+GMR
</td>
<td style="text-align:left;">
In-sample
</td>
<td style="text-align:right;">
82
</td>
<td style="text-align:right;">
82
</td>
<td style="text-align:right;">
57
</td>
<td style="text-align:right;">
53
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+GMR
</td>
<td style="text-align:left;">
Out-of-sample
</td>
<td style="text-align:right;">
71
</td>
<td style="text-align:right;">
71
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
12
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+DiCER
</td>
<td style="text-align:left;">
In-sample
</td>
<td style="text-align:right;">
82
</td>
<td style="text-align:right;">
82
</td>
<td style="text-align:right;">
57
</td>
<td style="text-align:right;">
53
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+DiCER
</td>
<td style="text-align:left;">
Out-of-sample
</td>
<td style="text-align:right;">
74
</td>
<td style="text-align:right;">
74
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
7
</td>
</tr>
</tbody>
</table>

This table summarises the number of ROIs for which raw accuracy or
balanced accuracy is significantly greater than the model-free shuffle
null distribution, both before and after adjusting for multiple
comparisons with BH-FDR.

### CV linear SVM – inv prob

![](Step3_catch22_README_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

I’ve plotted the distribution of null accuracies (teal) alongside the
actual accuracies (pink) for the 82 ROIs on the left. Let’s zoom in on
AROMA+2P and pick the five brain regions with the highest
cross-validated balanced accuracy:

![](Step3_catch22_README_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Noise_Proc
</th>
<th style="text-align:left;">
Sample_Type
</th>
<th style="text-align:right;">
num_sig_acc
</th>
<th style="text-align:right;">
num_sig_acc_fdr
</th>
<th style="text-align:right;">
num_sig_bacc
</th>
<th style="text-align:right;">
num_sig_bacc_fdr
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
AROMA+2P
</td>
<td style="text-align:left;">
In-sample
</td>
<td style="text-align:right;">
73
</td>
<td style="text-align:right;">
64
</td>
<td style="text-align:right;">
82
</td>
<td style="text-align:right;">
82
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P
</td>
<td style="text-align:left;">
Out-of-sample
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
18
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+GMR
</td>
<td style="text-align:left;">
In-sample
</td>
<td style="text-align:right;">
76
</td>
<td style="text-align:right;">
73
</td>
<td style="text-align:right;">
82
</td>
<td style="text-align:right;">
82
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+GMR
</td>
<td style="text-align:left;">
Out-of-sample
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
25
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+DiCER
</td>
<td style="text-align:left;">
In-sample
</td>
<td style="text-align:right;">
77
</td>
<td style="text-align:right;">
72
</td>
<td style="text-align:right;">
82
</td>
<td style="text-align:right;">
82
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+DiCER
</td>
<td style="text-align:left;">
Out-of-sample
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
24
</td>
</tr>
</tbody>
</table>

This table summarises the number of ROIs for which raw accuracy or
balanced accuracy is significantly greater than the model-free shuffle
null distribution, both before and after adjusting for multiple
comparisons with BH-FDR.

### CV linear SVM – SMOTE

![](Step3_catch22_README_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

I’ve plotted the distribution of null accuracies (teal) alongside the
actual accuracies (pink) for the 82 ROIs on the left. Let’s zoom in on
AROMA+2P and pick the five brain regions with the highest
cross-validated balanced accuracy:

![](Step3_catch22_README_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Noise_Proc
</th>
<th style="text-align:left;">
Sample_Type
</th>
<th style="text-align:right;">
num_sig_acc
</th>
<th style="text-align:right;">
num_sig_acc_fdr
</th>
<th style="text-align:right;">
num_sig_bacc
</th>
<th style="text-align:right;">
num_sig_bacc_fdr
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
AROMA+2P
</td>
<td style="text-align:left;">
In-sample
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
82
</td>
<td style="text-align:right;">
82
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P
</td>
<td style="text-align:left;">
Out-of-sample
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
16
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+GMR
</td>
<td style="text-align:left;">
In-sample
</td>
<td style="text-align:right;">
82
</td>
<td style="text-align:right;">
82
</td>
<td style="text-align:right;">
82
</td>
<td style="text-align:right;">
82
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+GMR
</td>
<td style="text-align:left;">
Out-of-sample
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
28
</td>
<td style="text-align:right;">
28
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+DiCER
</td>
<td style="text-align:left;">
In-sample
</td>
<td style="text-align:right;">
82
</td>
<td style="text-align:right;">
82
</td>
<td style="text-align:right;">
82
</td>
<td style="text-align:right;">
82
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+DiCER
</td>
<td style="text-align:left;">
Out-of-sample
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
21
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
over N=10 iterations per ROI (N=82) and pooling the resulting accuracy
and balanced accuracy values, to generate empirical null distributions
of N=820 data points each, respectively.

### Unweighted

![](Step3_catch22_README_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

The fitted empirical null model distribution is fairly similar to the
real accuracy and balanced accuracy values using in-sample linear SVM
with no reweighting.

![](Step3_catch22_README_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Noise_Proc
</th>
<th style="text-align:left;">
Sample_Type
</th>
<th style="text-align:right;">
num_sig_acc
</th>
<th style="text-align:right;">
num_sig_acc_fdr
</th>
<th style="text-align:right;">
num_sig_bacc
</th>
<th style="text-align:right;">
num_sig_bacc_fdr
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
AROMA+2P
</td>
<td style="text-align:left;">
In-sample
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P
</td>
<td style="text-align:left;">
Out-of-sample
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+GMR
</td>
<td style="text-align:left;">
In-sample
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+GMR
</td>
<td style="text-align:left;">
Out-of-sample
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+DiCER
</td>
<td style="text-align:left;">
In-sample
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+DiCER
</td>
<td style="text-align:left;">
Out-of-sample
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
3
</td>
</tr>
</tbody>
</table>

### Inverse probability weighted

![](Step3_catch22_README_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

The fitted empirical null model distribution is fairly similar to the
real accuracy and balanced accuracy values using in-sample linear SVM
with no reweighting.

![](Step3_catch22_README_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Noise_Proc
</th>
<th style="text-align:left;">
Sample_Type
</th>
<th style="text-align:right;">
num_sig_acc
</th>
<th style="text-align:right;">
num_sig_acc_fdr
</th>
<th style="text-align:right;">
num_sig_bacc
</th>
<th style="text-align:right;">
num_sig_bacc_fdr
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
AROMA+2P
</td>
<td style="text-align:left;">
In-sample
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P
</td>
<td style="text-align:left;">
Out-of-sample
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+GMR
</td>
<td style="text-align:left;">
In-sample
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
10
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+GMR
</td>
<td style="text-align:left;">
Out-of-sample
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+DiCER
</td>
<td style="text-align:left;">
In-sample
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+DiCER
</td>
<td style="text-align:left;">
Out-of-sample
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
2
</td>
</tr>
</tbody>
</table>

### SMOTE

![](Step3_catch22_README_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

The fitted empirical null model distribution is fairly similar to the
real accuracy and balanced accuracy values using in-sample linear SVM
with no reweighting.

![](Step3_catch22_README_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

<table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Noise_Proc
</th>
<th style="text-align:left;">
Sample_Type
</th>
<th style="text-align:right;">
num_sig_acc
</th>
<th style="text-align:right;">
num_sig_acc_fdr
</th>
<th style="text-align:right;">
num_sig_bacc
</th>
<th style="text-align:right;">
num_sig_bacc_fdr
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
AROMA+2P
</td>
<td style="text-align:left;">
In-sample
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P
</td>
<td style="text-align:left;">
Out-of-sample
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+GMR
</td>
<td style="text-align:left;">
In-sample
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+GMR
</td>
<td style="text-align:left;">
Out-of-sample
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+DiCER
</td>
<td style="text-align:left;">
In-sample
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
AROMA+2P+DiCER
</td>
<td style="text-align:left;">
Out-of-sample
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
2
</td>
</tr>
</tbody>
</table>

## Comparing model-free shuffle with pooled empirical null distributions

![](Step3_catch22_README_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->
