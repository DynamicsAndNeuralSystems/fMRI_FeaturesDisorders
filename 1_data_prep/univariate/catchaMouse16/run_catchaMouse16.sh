#!/usr/bin/env bash

study=/media/sf_Shared_Folder/PhD_work/data/scz/UCLA/catchaMouse16_data
catchaMouse_dir=/media/sf_Shared_Folder/github/catchaMouse16/C

cd $catchaMouse_dir

input_dir=$study/input_csv
output_dir=$study/output_csv

for input_file in $(dir -1 $input_dir | grep "_TS.csv"); do
  basename=$(basename $input_file _TS.csv)
  output_file=$output_dir/${basename}_cm16.csv

  ./run_feat $input_dir/$input_file $output_file
done
