##########################################################################################
cd $github_dir/fMRI_FeaturesDisorders/prep_data_and_QC/

# ##########################################################################################
# Read metadata into R and python
echo "Now preparing metadata"
Rscript dataset_specific_files/parse_metadata_${dataset_ID}.R

# ##########################################################################################
# Pre-filter subjects based on head movement using the lenient threshold from Parkes et al. (2018)
echo "Now finding subjects to drop based on head movement"
Rscript filter_subjects_based_on_head_movement.R --dataset_ID $dataset_ID \
--data_path $data_path \
--sample_metadata_file $sample_metadata_file

# ##########################################################################################
# Prep univariate data in R
echo "Now preparing univariate data"
univariate_feature_set=catch25
qsub -v dataset_ID=$dataset_ID,data_path=$data_path,univariate_feature_set=$univariate_feature_set,sample_metadata_file=$sample_metadata_file,brain_region_lookup=$brain_region_lookup,noise_proc=$noise_proc \
-N prepare_univariate_data_${dataset_ID} \
-o $github_dir/fMRI_FeaturesDisorders/cluster_output/prepare_univariate_data_${dataset_ID}_${feature_set}_out.txt \
-m a \
call_prepare_univariate_data.pbs

########################################################################################
# Prep pairwise data
# Get data into .npy files
echo "Now preparing pairwise data"
qsub -N prepare_pairwise_data_${dataset_ID} -q yossarian -j oe \
-v github_dir=$github_dir,conda_env=$conda_env,data_path=$data_path,noise_proc=$noise_proc,dataset_ID=$dataset_ID,sample_metadata_file=$sample_metadata_file,label_vars=$label_vars \
-o $github_dir/fMRI_FeaturesDisorders/cluster_output/prepare_pairwise_data_${dataset_ID}.txt \
-l select=1:ncpus=1:mem=20GB -l walltime=2:00:00 -M $email -m a -V \
prepare_pairwise_data.sh

# Run pyspi-distribute
echo "Now submitting pyspi-distribute jobs"
bash $github_dir/fMRI_FeaturesDisorders/prep_data_and_QC/call_run_pyspi_distribute.sh

