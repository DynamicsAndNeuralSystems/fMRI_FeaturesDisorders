# Feature extraction and classification for schizophrenia neuroimaging

This repository provides code that facilitates:

a) feature extraction using [_catch22_](https://github.com/chlubba/catch22) from a Matlab .mat file containing time-series data


### 1. Generating the time-series catch22 feature matrix

First run `prepare_TS_feature_data.R`

e.g. `Rscript prepare_TS_feature_data.R --mat_file D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/new/UCLA_time_series_four_groups.mat --label_metadata D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/UCLA/participants.csv --rdata_path D:/Virtual_Machines/Shared_Folder/PhD_work/data/scz/Rdata/`

Include the --overwrite flag if you want to overwrite the time series .rds objects if they already exist.


To calculate the total number of controls and patients (SCZ), use the below code:

`get_dx_breakdown(path_to_subject_csv_file)`


### 2. Analysis functions and their outputs

#### a) Region-by-Region Analysis

To look at the classification accuracy for e.g. the left entorhinal cortex, you could run:

`region_by_region_analysis(ROI="ctx-lh-entorhinal", feature_matrix = feature_matrix, display_figures = T, return_restable=F)`