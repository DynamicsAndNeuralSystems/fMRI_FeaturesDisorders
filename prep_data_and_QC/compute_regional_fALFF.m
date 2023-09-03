UCLA_CNP_data_path = '/Users/abry4213/data/UCLA_CNP/';
ABIDE_ASD_data_path = '/Users/abry4213/data/ABIDE_ASD/';

% Define TR as 2s
samplingPeriod = 2;

% Load in feather data file
UCLA_CNP_TS_file = "UCLA_CNP_AROMA_2P_GMR_fMRI_TS.feather";

% Path to the python script that converts the .feather file to .mat if
% needed
python_script_path = "/Users/abry4213/Library/CloudStorage/OneDrive-TheUniversityofSydney(Students)/github/fMRI_FeaturesDisorders/prep_data_and_QC/feather_to_mat.py";

% Command to run the python script from MATLAB
feather_file_base = "/Users/abry4213/data/UCLA_CNP/raw_data/UCLA_CNP_AROMA_2P_GMR_fMRI_TS";
command = sprintf('/Users/abry4213/anaconda3/envs/pyspi/bin/python3 "%s" "%s"', python_script_path, feather_file_base);
[~, result] = system(command);
result

% Load the .mat file
load(sprintf('%s.mat', feather_file_base));

[m, n] = size(feather_to_mat_struct.Sample_ID);
% Initialise fALFF results
fALFF_data = cell(m,3);

for i = 1:m
    sample_ID = feather_to_mat_struct.Sample_ID(i, :); % Extract sample ID
    brain_region = feather_to_mat_struct.Brain_Region(i, :); % Extract brain region
    time_series = transpose(feather_to_mat_struct.Time_Series(i, :)); % Extract time series

    % Compute fALFF for this participant/region
    SP_fALFF_res = SP_fALFF(time_series, samplingPeriod);
    
    fALFF = SP_fALFF_res.fALFF;
    ALFF = SP_fALFF_res.ALFF;

    % Add results to mat
    fALFF_data{i, 1} = sample_ID;
    fALFF_data{i, 2} = brain_region;
    fALFF_data{i, 3} = fALFF;
    fALFF_data{i, 4} = ALFF;
end


save(strcat(UCLA_CNP_data_path, 'processed_data/UCLA_CNP_AROMA_2P_GMR_ALFF_fALFF.mat'), 'fALFF_data', '-v6')