
# UCLA Schizophrenia
Rscript driver_script_data_prep.R --project_path /project/hctsa/annie/github/ \
--github_dir /project/hctsa/annie/github/fMRI_FeaturesDisorders/ \
--data_path /project/hctsa/annie/data/UCLA_Schizophrenia/ \
--rdata_path /project/hctsa/annie/data/UCLA_Schizophrenia/Rdata/ \
--univariate_feature_set catch22 --parcellation_name aparc+aseg \
--input_mat_file new/UCLA_time_series_four_groups.mat \
--subject_csv participants.csv \
--brain_region_lookup \
--noise_procs AROMA+2P AROMA+2P+GMR AROMA+2P+DICER \
--dataset_ID UCLA_Schizophrenia 

# ABIDE ASD
Rscript driver_script_data_prep.R --project_path /project/hctsa/annie/github/ \
--github_dir /project/hctsa/annie/github/fMRI_FeaturesDisorders/ \
--data_path /project/hctsa/annie/data/ABIDE_ASD/ \
--rdata_path /project/hctsa/annie/data/ABIDE_ASD/Rdata/ \
--univariate_feature_set catch22 \
--parcellation_name harvard_oxford_cort_prob_2mm \
--input_mat_file NA \
--subject_csv participants.csv \
--brain_region_lookup Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv \
--noise_procs FC1000 \
--dataset_ID ABIDE_ASD 