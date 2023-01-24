
##########################################################################################
cd $github_dir/fMRI_FeaturesDisorders/prep_data_and_QC/

# #########################################################################################
# Clean up pairwise data
echo "Now preparing pairwise data"
for univariate_feature_set in catch2 catch22 catch24; do 
    qsub -N clean_pairwise_data_${dataset_ID}_${pairwise_feature_set} -q yossarian -j oe \
    -v github_dir=$github_dir,conda_env=$conda_env,data_path=$data_path,noise_proc=$noise_proc,dataset_ID=$dataset_ID,pairwise_feature_set=$pairwise_feature_set,univariate_feature_set=$univariate_feature_set,brain_region_lookup=$brain_region_lookup \
    -o $github_dir/fMRI_FeaturesDisorders/cluster_output/clean_pairwise_data_${dataset_ID}_${pairwise_feature_set}.txt \
    -l select=1:ncpus=1:mem=20GB -l walltime=1:00:00 -M $email -m a -V \
    call_clean_pairwise_data.pbs
done