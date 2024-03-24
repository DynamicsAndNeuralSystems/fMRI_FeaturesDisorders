#!/usr/bin/env bash

parsing_python_script=/headnode1/abry4213/github/PISA_analysis/classification_analysis/parse_npy_for_classification.py
data_path=/headnode1/abry4213/data/TS_feature_manuscript
amyloid_thresh_type=full
output_data_path=/headnode1/abry4213/data/TS_feature_manuscript/classification_results
num_folds=10
num_repeats=10
mem=10
queue=defaultQ
classifier_type=Linear_SVM_sklearn
num_jobs=4
num_null_iters=200

# Specific here
dataset_ID=ABIDE
disorder=ASD
analysis_type=Brain_Region
class_labels_file=${data_path}/input_data/${dataset_ID}_${disorder}_dx_labels.npy
sample_IDs_file=${data_path}/input_data/${dataset_ID}_${disorder}_sample_IDs.npy

core_count=0

cat ${data_path}/time_series_features/processed_numpy_files/${dataset_ID}_${disorder}_models.txt | while read line || [[ -n $line ]]
do
    loop_count=0
    model_name=$line 
    input_npy_file=/headnode1/abry4213/data/TS_feature_manuscript/time_series_features/processed_numpy_files/${model_name}.npy

    for null_start in 1 201 401 601 801; do 
        ((loop_count+=1))

        # Check if nulls were already done
        if test -f /headnode1/abry4213/data/TS_feature_manuscript/classification_results/null_results/${model_name}_${classifier_type}_${num_repeats}_repeats_${num_folds}_folds_CV_nulls_from_${null_start}.feather; then
            continue
        fi

        sleep 0.3
        ((core_count+=2))
        # Check if core_count has reached 30
        if [ $core_count -ge 81 ]; then
            core_count=0
            sleep 30m
        fi

        echo $model_name $null_start; 
        # Make the cluster job results output directory
        mkdir -p /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/${dataset_ID}_${disorder}_classification/

        # Submit the classification script
        cmd="qsub -o /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/${dataset_ID}_${disorder}_classification/${model_name}_${classifier_type}_${null_start}.out \
                -N ${model_name}_${null_start}_classification \
                -q ${queue} \
                -l select=1:ncpus=${num_jobs}:mem=${mem}GB:mpiprocs=${num_jobs} \
                -v dataset_ID=$dataset_ID,disorder=$disorder,input_npy_file=$input_npy_file,class_labels_file=$class_labels_file,sample_IDs_file=$sample_IDs_file,output_data_path=$output_data_path,analysis_type=$analysis_type,feature_name=$model_name,classifier_type=$classifier_type,num_folds=$num_folds,num_repeats=$num_repeats,num_jobs=$num_jobs,num_null_iters=$num_null_iters,null_starting_point=$null_start,random_num_for_perm=$loop_count \
                run_classification_for_model.pbs"
        $cmd
    done
done

# Merge univariate null results
python3 merge_null_results.py --dataset_ID ${dataset_ID} --disorder ${disorder} --data_path ${output_data_path}/null_results --main_analysis_type Univariate_catch25 --classifier_type ${classifier_type}