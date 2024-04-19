#!/usr/bin/env bash
num_running_cores=3
num_cores=1
num_jobs=2
mem_univariate=2
mem_undirected=7
num_jobs_directed=6
mem_directed=15
data_path=/headnode1/abry4213/data/TS_feature_manuscript
queue=yossarian
num_null_iters=1000

#################################################################################
# Univariate analyses
num_cores=1
num_jobs=2
num_GB=5
dataset_ID=UCLA_CNP
queue=defaultQ

for disorder in SCZ BP ADHD; do
    cmd="qsub -o /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/all_univariate_classification_${disorder}_^array_index^.out \
        -N ${disorder}_^array_index^_catch25 \
        -l select=1:ncpus=${num_cores}:mem=${num_GB}GB:mpiprocs=${num_cores} \
        -v num_null_iters=${num_null_iters},num_jobs=${num_jobs},dataset_ID=${dataset_ID},disorder=${disorder},input_model_file=${data_path}/time_series_features/processed_numpy_files/${dataset_ID}_${disorder}_univariate_models.txt \
        -J 1-36 \
        -q $queue \
        array_for_classification.pbs"
    $cmd
    sleep 3m
    cmd="qsub -o /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/all_univariate_classification_${disorder}_^array_index^.out \
        -N ${disorder}_^array_index^_catch25 \
        -l select=1:ncpus=${num_cores}:mem=${num_GB}GB:mpiprocs=${num_cores} \
        -v num_null_iters=${num_null_iters},num_jobs=${num_jobs},dataset_ID=${dataset_ID},disorder=${disorder},input_model_file=${data_path}/time_series_features/processed_numpy_files/${dataset_ID}_${disorder}_univariate_models.txt \
        -J 37-72 \
        -q $queue \
        array_for_classification.pbs"
    $cmd
    sleep 3m
    cmd="qsub -o /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/all_univariate_classification_${disorder}_^array_index^.out \
        -N ${disorder}_^array_index^_catch25 \
        -l select=1:ncpus=${num_cores}:mem=${num_GB}GB:mpiprocs=${num_cores} \
        -v num_null_iters=${num_null_iters},num_jobs=${num_jobs},dataset_ID=${dataset_ID},disorder=${disorder},input_model_file=${data_path}/time_series_features/processed_numpy_files/${dataset_ID}_${disorder}_univariate_models.txt \
        -J 73-107 \
        -q $queue \
        array_for_classification.pbs"
    $cmd
    cmd="qsub -o /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/all_univariate_classification_${disorder}_^array_index^.out \
        -N ${disorder}_^array_index^_catch25 \
        -l select=1:ncpus=${num_cores}:mem=${num_GB}GB:mpiprocs=${num_cores} \
        -v line_to_read=108,num_null_iters=${num_null_iters},num_jobs=${num_jobs},dataset_ID=${dataset_ID},disorder=${disorder},input_model_file=${data_path}/time_series_features/processed_numpy_files/${dataset_ID}_${disorder}_univariate_models.txt \
        -q $queue \
        array_for_classification.pbs"
    $cmd
    sleep 3m
done

for disorder in ASD; do
    dataset_ID=ABIDE
    cmd="qsub -o /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/all_univariate_classification_${disorder}_^array_index^.out \
        -N ${disorder}_catch25 \
        -l select=1:ncpus=${num_cores}:mem=${num_GB}GB:mpiprocs=${num_cores} \
        -v num_null_iters=${num_null_iters},num_jobs=${num_jobs},dataset_ID=${dataset_ID},disorder=${disorder},input_model_file=${data_path}/time_series_features/processed_numpy_files/${dataset_ID}_${disorder}_univariate_models.txt \
        -J 1-37 \
        -q $queue \
        array_for_classification.pbs"
    $cmd
    sleep 5m
    cmd="qsub -o /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/all_univariate_classification_${disorder}_^array_index^.out \
        -N ${disorder}_catch25 \
        -l select=1:ncpus=${num_cores}:mem=${num_GB}GB:mpiprocs=${num_cores} \
        -v num_null_iters=${num_null_iters},num_jobs=${num_jobs},dataset_ID=${dataset_ID},disorder=${disorder},input_model_file=${data_path}/time_series_features/processed_numpy_files/${dataset_ID}_${disorder}_univariate_models.txt \
        -J 38-73 \
        -q $queue \
        array_for_classification.pbs"
    $cmd
    cmd="qsub -o /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/all_univariate_classification_${disorder}_^array_index^.out \
        -N ${disorder}_catch25 \
        -l select=1:ncpus=${num_cores}:mem=${num_GB}GB:mpiprocs=${num_cores} \
        -v line_to_read=74,num_null_iters=${num_null_iters},num_jobs=${num_jobs},dataset_ID=${dataset_ID},disorder=${disorder},input_model_file=${data_path}/time_series_features/processed_numpy_files/${dataset_ID}_${disorder}_univariate_models.txt \
        -q $queue \
        array_for_classification.pbs"
    $cmd
done

# # Merge results
# for disorder in SCZ BP ADHD; do
#     python3 merge_main_and_null_results.py --dataset_ID UCLA_CNP --disorder $disorder \
#             --main_analysis_type Univariate_catch25 --classifier_type Linear_SVM_sklearn \
#             --data_path /headnode1/abry4213/data/TS_feature_manuscript/classification_results
# done
# for disorder in ASD; do
#     python3 merge_main_and_null_results.py --dataset_ID ABIDE --disorder $disorder \
#             --main_analysis_type Univariate_catch25 --classifier_type Linear_SVM_sklearn \
#             --data_path /headnode1/abry4213/data/TS_feature_manuscript/classification_results
# done

#################################################################################

# Pairwise
num_running_cores=3
num_cores=1
num_jobs=2
mem_univariate=2
mem_undirected=7
num_jobs_directed=6
mem_directed=15
queue=defaultQ
num_null_iters=1000

# for disorder in SCZ BP ADHD; do
    # # 1-3
    # queue=defaultQ
    # cmd="qsub -o /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/pyspi14_SPI_classification_${disorder}.out \
    #     -N ${disorder}_pyspi14 \
    #     -l select=1:ncpus=${num_cores}:mem=${mem_undirected}GB:mpiprocs=${num_cores} \
    #     -v num_null_iters=${num_null_iters},num_jobs=${num_jobs},dataset_ID=${dataset_ID},disorder=${disorder},input_model_file=${data_path}/time_series_features/processed_numpy_files/${dataset_ID}_${disorder}_pairwise_models.txt \
    #     -J 1-3 \
    #     -q $queue \
    #     array_for_classification.pbs"
    # $cmd

    # # 7-8
    # queue=defaultQ
    # cmd="qsub -o /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/pyspi14_SPI_classification_${disorder}.out \
    #     -N ${disorder}_pyspi14 \
    #     -l select=1:ncpus=${num_cores}:mem=${mem_undirected}GB:mpiprocs=${num_cores} \
    #     -v num_null_iters=${num_null_iters},num_jobs=${num_jobs},dataset_ID=${dataset_ID},disorder=${disorder},input_model_file=${data_path}/time_series_features/processed_numpy_files/${dataset_ID}_${disorder}_pairwise_models.txt \
    #     -J 7-8 \
    #     -q $queue \
    #     array_for_classification.pbs"
    # $cmd

    # # 12-14
    # queue=yossarian
    # cmd="qsub -o /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/pyspi14_SPI_classification_${disorder}.out \
    #     -N ${disorder}_pyspi14 \
    #     -l select=1:ncpus=${num_cores}:mem=${mem_undirected}GB:mpiprocs=${num_cores} \
    #     -v num_null_iters=${num_null_iters},num_jobs=${num_jobs},dataset_ID=${dataset_ID},disorder=${disorder},input_model_file=${data_path}/time_series_features/processed_numpy_files/${dataset_ID}_${disorder}_pairwise_models.txt \
    #     -J 12-14 \
    #     -q $queue \
    #     array_for_classification.pbs"
    # $cmd

    # # 4-6
    # queue=yossarian
    # cmd="qsub -o /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/pyspi14_SPI_classification_${disorder}.out \
    #     -N ${disorder}_pyspi14 \
    #     -l select=1:ncpus=${num_jobs_directed}:mem=${mem_directed}GB:mpiprocs=${num_jobs_directed} \
    #     -v num_null_iters=${num_null_iters},num_jobs=${num_jobs_directed},dataset_ID=${dataset_ID},disorder=${disorder},input_model_file=${data_path}/time_series_features/processed_numpy_files/${dataset_ID}_${disorder}_pairwise_models.txt \
    #     -J 4-6 \
    #     -q $queue \
    #     array_for_classification.pbs"
    # $cmd

#     # 9-11
#     cmd="qsub -o /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/pyspi14_SPI_classification_${disorder}.out \
#         -N ${disorder}_pyspi14 \
#         -l select=1:ncpus=${num_jobs_directed}:mem=${mem_directed}GB:mpiprocs=${num_jobs_directed} \
#         -v num_null_iters=${num_null_iters},num_jobs=${num_jobs_directed},dataset_ID=${dataset_ID},disorder=${disorder},input_model_file=${data_path}/time_series_features/processed_numpy_files/${dataset_ID}_${disorder}_pairwise_models.txt \
#         -J 9-11 \
#         -q $queue \
#         array_for_classification.pbs"
#     $cmd
# done

# dataset_ID=ABIDE
# num_jobs=6
# num_cores=30
# num_null_iters=1000
# queue=yossarian
# for disorder in ASD; do
#     cmd="qsub -o /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/pyspi14_SPI_classification_${disorder}.out \
#         -N ${disorder}_pyspi14 \
#         -l select=1:ncpus=${num_jobs}:mem=${num_cores}GB:mpiprocs=${num_jobs} \
#         -v num_null_iters=${num_null_iters},num_jobs=${num_jobs},dataset_ID=${dataset_ID},disorder=${disorder},input_model_file=${data_path}/time_series_features/processed_numpy_files/${dataset_ID}_${disorder}_pairwise_models.txt \
#         -J 1-14 \
#         -q $queue \
#         array_for_classification.pbs"
#     $cmd
# done

#################################################################################
# Combined univariate+pairwise
# num_cores=6
# mem_directed=15
# num_null_iters=1000
# queue=yossarian

# dataset_ID=UCLA_CNP
# for disorder in SCZ BP ADHD; do
#     cmd="qsub -o /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/combined_classification_${disorder}.out \
#         -N ${disorder}_pyspi14 \
#         -l select=1:ncpus=${num_cores}:mem=${mem_directed}GB:mpiprocs=${num_cores} \
#         -v num_null_iters=${num_null_iters},num_jobs=${num_cores},dataset_ID=${dataset_ID},disorder=${disorder},input_model_file=${data_path}/time_series_features/processed_numpy_files/${dataset_ID}_${disorder}_combined_univariate_pairwise_models.txt \
#         -J 1-14 \
#         -q $queue \
#         array_for_classification.pbs"
#     $cmd
#     sleep 10h
# done

# dataset_ID=ABIDE
# queue=yossarian
# num_cores=6
# mem_GB=15
# num_null_iters=1000

# for disorder in ASD; do
#     cmd="qsub -o /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/combined_classification_${disorder}.out \
#         -N ${disorder}_pyspi14 \
#         -l select=1:ncpus=${num_cores}:mem=${mem_GB}GB:mpiprocs=${num_cores} \
#         -v num_null_iters=${num_null_iters},num_jobs=${num_cores},dataset_ID=${dataset_ID},disorder=${disorder},input_model_file=${data_path}/time_series_features/processed_numpy_files/${dataset_ID}_${disorder}_combined_univariate_pairwise_models.txt \
#         -J 1-14 \
#         -q $queue \
#         array_for_classification.pbs"
#     $cmd
# done

