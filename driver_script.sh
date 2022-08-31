##########################################################################################
# Prep univariate data
qsub -v run_number=1 \
-N prepare_univariate_data1 \
-o /headnode/abry4213/github/fMRI_FeaturesDisorders/cluster_output/prepare_univariate_data1_out.txt \
data_prep_and_QC/call_prepare_univariate_data.pbs

qsub -v run_number=2 \
-N prepare_univariate_data2 \
-o /headnode/abry4213/github/fMRI_FeaturesDisorders/cluster_output/prepare_univariate_data2_out.txt \
data_prep_and_QC/call_prepare_univariate_data.pbs

qsub -v run_number=3 \
-N prepare_univariate_data3 \
-o /headnode/abry4213/github/fMRI_FeaturesDisorders/cluster_output/prepare_univariate_data3_out.txt \
data_prep_and_QC/call_prepare_univariate_data.pbs

qsub -v run_number=4 \
-N prepare_univariate_data4 \
-o /headnode/abry4213/github/fMRI_FeaturesDisorders/cluster_output/prepare_univariate_data4_out.txt \
data_prep_and_QC/call_prepare_univariate_data.pbs

qsub -v run_number=5 \
-N prepare_univariate_data5 \
-o /headnode/abry4213/github/fMRI_FeaturesDisorders/cluster_output/prepare_univariate_data5_out.txt \
data_prep_and_QC/call_prepare_univariate_data.pbs

# Prep pairwise data
# Get data into .npy files
qsub data_prep_and_QC/call_prepare_pairwise_data.pbs 
# Run pyspi-distribute
bash data_prep_and_QC/call_run_pyspi_distribute.sh 
# Integrate results from pyspi-distribute
qsub data_prep_and_QC/call_clean_pairwise_data.pbs 

# Merge subjects with univariate + pairwise data
qsub data_prep_and_QC/call_merge_samples_univariate_pairwise.pbs

##########################################################################################
# Univariate linear SVM
qsub univariate_analysis/call_univariate_classification.pbs 

# Generate null model fits
bash univariate_analysis/call_univariate_null_model_generation.sh

# Integrate null model fits and calculate p-values
qsub univariate_analysis/call_univariate_null_model_analysis.pbs

##########################################################################################
# Pairwise linear SVM
qsub pairwise_analysis/call_pairwise_classification.pbs 
