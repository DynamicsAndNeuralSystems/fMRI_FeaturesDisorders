##########################################################################################
# Prep univariate data
qsub data_prep_and_QC/call_prepare_univariate_data.pbs

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
