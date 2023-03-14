##########################################################################################
export github_dir=/headnode1/abry4213/github/
export pairwise_feature_set="pyspi14"
export univariate_feature_set="catch22"
export email="abry4213@uni.sydney.edu.au"
export conda_env="pyspi"
export cluster_queue="yossarian"
export python_to_use=/headnode1/abry4213/.conda/envs/${conda_env}/bin/python3
export SPI_directionality_file=${github_dir}/fMRI_FeaturesDisorders/classification_analysis/SPI_Direction_Info.csv
export scaling_type="mixedsigmoid"
export num_folds=10
export num_repeats=10
export num_null_iters=1000

#module load Anaconda3-5.1.0
/usr/physics/Modules/3.2.8/bin/modulecmd bash load Anaconda3-5.1.0 --silent
source /usr/physics/python/anaconda3/etc/profile.d/conda.sh 

# Activate the given conda environment
conda activate $conda_env

##########################################################################################
# UCLA CNP

# Define data paths
export dataset_ID="UCLA_CNP"
export data_path=/headnode1/abry4213/data/UCLA_CNP/
export sample_metadata_file="UCLA_CNP_sample_metadata.feather"
export noise_proc="AROMA+2P+GMR"
export noise_label="AROMA_2P_GMR"
export brain_region_lookup="UCLA_CNP_Brain_Region_Lookup.feather"
# Run univariate and pairwise linear SVM

for comparison_group in Schizophrenia Bipolar ADHD; do
    pairwise_feature_file=$data_path/processed_data/${dataset_ID}_${noise_label}_${pairwise_feature_set}_filtered.feather
    univariate_feature_file=$data_path/processed_data/${dataset_ID}_${noise_label}_${univariate_feature_set}_filtered.feather

    # Run univariate SVM

    # Main SVM
    num_jobs=10
    job_memory=20
    num_hours=3
    qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,univariate_feature_file=$univariate_feature_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_folds=$num_folds,num_jobs=$num_jobs,num_repeats=$num_repeats \
    -N ${dataset_ID}_${comparison_group}_${univariate_feature_set} \
    -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_univariate_classification_${dataset_ID}_${comparison_group}_${univariate_feature_set}_${scaling_type}_scaler_out.txt \
    -m a -M $email -l select=1:ncpus=$num_jobs:mem=${job_memory}GB,walltime=${num_hours}:00:00 -q $cluster_queue \
    $github_dir/fMRI_FeaturesDisorders/classification_analysis/call_univariate_classification.pbs

    # # Nulls
    # job_memory=20
    # num_hours=24
    # num_jobs=10
    # qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,univariate_feature_file=$univariate_feature_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_null_iters=$num_null_iters,num_folds=$num_folds,num_jobs=$num_jobs,num_repeats=$num_repeats \
    # -N ${dataset_ID}_${comparison_group}_${univariate_feature_set} \
    # -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_univariate_nulls_${dataset_ID}_${comparison_group}_${univariate_feature_set}_${scaling_type}_scaler_out.txt \
    # -m a -M $email -l select=1:ncpus=$num_jobs:mem=${job_memory}GB,walltime=${num_hours}:00:00 -q $cluster_queue \
    # $github_dir/fMRI_FeaturesDisorders/classification_analysis/call_univariate_nulls.pbs

    # Run movement-based univariate SVM
    job_memory=20
    num_jobs=10
    num_hours=3
    num_folds=5
    qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,univariate_feature_file=$univariate_feature_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_null_iters=$num_null_iters,num_folds=$num_folds,num_jobs=$num_jobs,num_repeats=$num_repeats \
    -N ${dataset_ID}_${comparison_group}_movement_${univariate_feature_set} \
    -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_movement_region_wise_SVM_${dataset_ID}_${comparison_group}_${univariate_feature_set}_${scaling_type}_scaler_out.txt \
    -m a -M $email -l select=1:ncpus=$num_jobs:mem=${job_memory}GB,walltime=${num_hours}:00:00 -q $cluster_queue \
    $github_dir/fMRI_FeaturesDisorders/prep_data_and_QC/movement/call_SVMs_by_movement.pbs

    # Pairwise Main SVM
    num_jobs=10
    job_memory=40
    num_hours=6
    num_folds=10
    qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,pairwise_feature_file=$pairwise_feature_file,SPI_directionality_file=$SPI_directionality_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_folds=$num_folds,num_jobs=$num_jobs,num_repeats=$num_repeats \
    -N ${dataset_ID}_${comparison_group}_${pairwise_feature_set} \
    -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_pairwise_classification_${dataset_ID}_${comparison_group}_${pairwise_feature_set}_${scaling_type}_scaler_out.txt \
    -m a -M $email -l select=1:ncpus=$num_jobs:mem=${job_memory}GB,walltime=${num_hours}:00:00 -q $cluster_queue \
    $github_dir/fMRI_FeaturesDisorders/classification_analysis/call_pairwise_classification.pbs

    # # Pairwise Nulls
    # num_jobs=10
    # job_memory=180
    # num_hours=96
    # qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,pairwise_feature_file=$pairwise_feature_file,SPI_directionality_file=$SPI_directionality_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_null_iters=$num_null_iters,num_folds=$num_folds,num_jobs=$num_jobs,num_repeats=$num_repeats \
    # -N ${dataset_ID}_${comparison_group}_${pairwise_feature_set} \
    # -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_pairwise_nulls_${dataset_ID}_${comparison_group}_${pairwise_feature_set}_${scaling_type}_scaler_out.txt \
    # -m a -M $email -l select=1:ncpus=$num_jobs:mem=${job_memory}GB,walltime=${num_hours}:00:00 -q $cluster_queue \
    # $github_dir/fMRI_FeaturesDisorders/classification_analysis/call_pairwise_nulls.pbs
    
    # Combined univariate + pairwise Main SVM
    num_jobs=10
    job_memory=80
    num_hours=4
    num_folds=10
    qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,univariate_feature_file=$univariate_feature_file,pairwise_feature_set=$pairwise_feature_set,pairwise_feature_file=$pairwise_feature_file,SPI_directionality_file=$SPI_directionality_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_jobs=$num_jobs,num_repeats=$num_repeats \
    -N ${dataset_ID}_${comparison_group}_${univariate_feature_set}_${pairwise_feature_set} \
    -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_combo_classification_${dataset_ID}_${comparison_group}_${univariate_feature_set}_${pairwise_feature_set}_${scaling_type}_scaler_out.txt \
    -m a -M $email -l select=1:ncpus=$num_jobs:mem=${job_memory}GB,walltime=${num_hours}:00:00 -q $cluster_queue  \
    $github_dir/fMRI_FeaturesDisorders/classification_analysis/call_combined_univariate_pairwise_classification.pbs

    # # Combined univariate + pairwise Nulls
    # num_jobs=10
    # job_memory=80
    # num_hours=120
    # qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,univariate_feature_file=$univariate_feature_file,pairwise_feature_set=$pairwise_feature_set,pairwise_feature_file=$pairwise_feature_file,SPI_directionality_file=$SPI_directionality_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_null_iters=$num_null_iter,num_jobs=$num_jobs,num_repeats=$num_repeats \
    # -N ${dataset_ID}_${comparison_group}_${univariate_feature_set}_${pairwise_feature_set} \
    # -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_combo_nulls_${dataset_ID}_${comparison_group}_${univariate_feature_set}_${pairwise_feature_set}_${scaling_type}_scaler_out.txt \
    # -m a -M $email -l select=1:ncpus=$num_jobs:mem=${job_memory}GB,walltime=${num_hours}:00:00 -q $cluster_queue  \
    # $github_dir/fMRI_FeaturesDisorders/classification_analysis/call_combined_univariate_pairwise_nulls.pbs
done

##########################################################################################
# ABIDE ASD

# Define data paths
export dataset_ID="ABIDE_ASD"
export data_path=/headnode1/abry4213/data/ABIDE_ASD/
export sample_metadata_file=ABIDE_ASD_sample_metadata.feather
export noise_proc="FC1000"
export noise_label="FC1000"
export brain_region_lookup="ABIDE_ASD_Brain_Region_Lookup.feather"

# Run univariate and pairwise linear SVM
for comparison_group in ASD; do
    univariate_feature_file=$data_path/processed_data/${dataset_ID}_${noise_label}_${univariate_feature_set}_filtered.feather
    pairwise_feature_file=$data_path/processed_data/${dataset_ID}_${noise_label}_${pairwise_feature_set}_filtered.feather
    # Run univariate SVM

    # Univariate Main SVM
    num_jobs=10
    job_memory=80
    num_hours=3
    qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,univariate_feature_file=$univariate_feature_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_folds=$num_folds,num_jobs=$num_jobs,num_repeats=$num_repeats \
    -N ${dataset_ID}_${comparison_group}_${univariate_feature_set} \
    -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_univariate_classification_${dataset_ID}_${comparison_group}_${univariate_feature_set}_${scaling_type}_scaler_out.txt \
    -m a -M $email -l select=1:ncpus=$num_jobs:mem=${job_memory}GB,walltime=${num_hours}:00:00 -q $cluster_queue \
    $github_dir/fMRI_FeaturesDisorders/classification_analysis/call_univariate_classification.pbs

    # # Univariate Nulls
    # num_jobs=10
    # job_memory=180
    # num_hours=48
    # qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,univariate_feature_file=$univariate_feature_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_null_iters=$num_null_iters,num_folds=$num_folds,num_jobs=$num_jobs,num_repeats=$num_repeats \
    # -N ${dataset_ID}_${comparison_group}_${univariate_feature_set} \
    # -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_univariate_nulls_${dataset_ID}_${comparison_group}_${univariate_feature_set}_${scaling_type}_scaler_out.txt \
    # -m a -M $email -l select=1:ncpus=$num_jobs:mem=${job_memory}GB,walltime=${num_hours}:00:00 -q $cluster_queue \
    # $github_dir/fMRI_FeaturesDisorders/classification_analysis/call_univariate_nulls.pbs

    # Run movement-based univariate SVM
    job_memory=80
    num_jobs=10
    num_hours=3
    num_folds=5
    qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,univariate_feature_file=$univariate_feature_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_null_iters=$num_null_iters,num_folds=$num_folds,num_jobs=$num_jobs,num_repeats=$num_repeats \
    -N ${dataset_ID}_${comparison_group}_movement_${univariate_feature_set} \
    -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_movement_region_wise_SVM_${dataset_ID}_${comparison_group}_${univariate_feature_set}_${scaling_type}_scaler_out.txt \
    -m a -M $email -l select=1:ncpus=$num_jobs:mem=${job_memory}GB,walltime=${num_hours}:00:00 -q $cluster_queue \
    $github_dir/fMRI_FeaturesDisorders/prep_data_and_QC/movement/call_SVMs_by_movement.pbs
        
    # Pairwise Main SVM
    num_jobs=10
    job_memory=80
    num_hours=6
    num_folds=10
    qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,pairwise_feature_file=$pairwise_feature_file,SPI_directionality_file=$SPI_directionality_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_folds=$num_folds,num_jobs=$num_jobs,num_repeats=$num_repeats \
    -N ${dataset_ID}_${comparison_group}_${pairwise_feature_set} \
    -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_pairwise_classification_${dataset_ID}_${comparison_group}_${pairwise_feature_set}_${scaling_type}_scaler_out.txt \
    -m a -M $email -l select=1:ncpus=$num_jobs:mem=${job_memory}GB,walltime=${num_hours}:00:00 -q $cluster_queue \
    $github_dir/fMRI_FeaturesDisorders/classification_analysis/call_pairwise_classification.pbs

    # # Pairwise Nulls
    # num_jobs=10
    # job_memory=180
    # num_hours=120
    # qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,pairwise_feature_file=$pairwise_feature_file,SPI_directionality_file=$SPI_directionality_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_null_iters=$num_null_iters,num_folds=$num_folds,num_jobs=$num_jobs,num_repeats=$num_repeats \
    # -N ${dataset_ID}_${comparison_group}_${pairwise_feature_set} \
    # -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_pairwise_nulls_${dataset_ID}_${comparison_group}_${pairwise_feature_set}_${scaling_type}_scaler_out.txt \
    # -m a -M $email -l select=1:ncpus=$num_jobs:mem=${job_memory}GB,walltime=${num_hours}:00:00 -q $cluster_queue \
    # $github_dir/fMRI_FeaturesDisorders/classification_analysis/call_pairwise_nulls.pbs

    # Combined univariate + pairwise Main SVM
    num_jobs=10
    job_memory=180
    num_hours=6
    num_folds=10
    qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,univariate_feature_file=$univariate_feature_file,pairwise_feature_set=$pairwise_feature_set,pairwise_feature_file=$pairwise_feature_file,SPI_directionality_file=$SPI_directionality_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_jobs=$num_jobs,num_repeats=$num_repeats \
    -N ${dataset_ID}_${comparison_group}_${univariate_feature_set}_${pairwise_feature_set} \
    -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_combo_classification_${dataset_ID}_${comparison_group}_${univariate_feature_set}_${pairwise_feature_set}_${scaling_type}_scaler_out.txt \
    -m a -M $email -l select=1:ncpus=$num_jobs:mem=${job_memory}GB,walltime=${num_hours}:00:00 -q $cluster_queue  \
    $github_dir/fMRI_FeaturesDisorders/classification_analysis/call_combined_univariate_pairwise_classification.pbs

    # # Combined univariate + pairwise Nulls
    # num_jobs=10
    # job_memory=180
    # num_hours=150
    # qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,univariate_feature_file=$univariate_feature_file,pairwise_feature_set=$pairwise_feature_set,pairwise_feature_file=$pairwise_feature_file,SPI_directionality_file=$SPI_directionality_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_null_iters=$num_null_iter,num_jobs=$num_jobs,num_repeats=$num_repeats \
    # -N ${dataset_ID}_${comparison_group}_${univariate_feature_set}_${pairwise_feature_set} \
    # -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_combo_nulls_${dataset_ID}_${comparison_group}_${univariate_feature_set}_${pairwise_feature_set}_${scaling_type}_scaler_out.txt \
    # -m a -M $email -l select=1:ncpus=$num_jobs:mem=${job_memory}GB,walltime=${num_hours}:00:00 -q $cluster_queue  \
    # $github_dir/fMRI_FeaturesDisorders/classification_analysis/call_combined_univariate_pairwise_nulls.pbs
done
