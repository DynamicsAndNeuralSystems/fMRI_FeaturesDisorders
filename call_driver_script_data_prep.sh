
# UCLA Schizophrenia
# Rscript driver_script_data_prep.R --project_path D:/Virtual_Machines/Shared_Folder/github/ \
# --github_dir D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/ \
# --data_path D:/Virtual_Machines/Shared_Folder/PhD_work/data/UCLA_Schizophrenia/ \
# --rdata_path D:/Virtual_Machines/Shared_Folder/PhD_work/data/UCLA_Schizophrenia/Rdata/ \
# --univariate_feature_set catch22 --parcellation_name aparc+aseg \
# --input_mat_file new/UCLA_time_series_four_groups.mat \
# --subject_csv participants.csv \
# --brain_region_lookup \
# --noise_procs AROMA+2P AROMA+2P+GMR AROMA+2P+DICER \
# --dataset_ID UCLA_Schizophrenia 

# ABIDE ASD
Rscript driver_script_data_prep.R --project_path D:/Virtual_Machines/Shared_Folder/github/ \
--github_dir D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/ \
--data_path D:/Virtual_Machines/Shared_Folder/PhD_work/data/ABIDE_ASD/ \
--rdata_path D:/Virtual_Machines/Shared_Folder/PhD_work/data/ABIDE_ASD/Rdata/ \
--univariate_feature_set catch22 \
--parcellation_name harvard_oxford_cort_prob_2mm \
--input_mat_file NA \
--subject_csv participants.csv \
--brain_region_lookup Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv \
--noise_procs FC1000 \
--dataset_ID ABIDE_ASD 