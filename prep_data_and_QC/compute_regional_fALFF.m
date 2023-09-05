UCLA_CNP_data_path = '/headnode1/abry4213/data/UCLA_CNP/';
ABIDE_ASD_data_path = '/headnode1/abry4213/data/ABIDE_ASD/';

% Add SP_fALFF script to path
addpath('/headnode1/abry4213/github/fMRI_FeaturesDisorders/prep_data_and_QC/');

% Define TR as 2s
samplingPeriod = 2;

% Load in feather data file
UCLA_CNP_TS_file = strcat(UCLA_CNP_data_path, "raw_data/UCLA_CNP_AROMA_2P_GMR_fMRI_TS.mat");

% Load the .mat file
load(UCLA_CNP_TS_file);

[m, n] = size(feather_to_mat_struct.Sample_ID);
% Initialise fALFF results
UCLA_CNP_fALFF_data = cell(m,3);

for i = 1:m
    sample_ID = feather_to_mat_struct.Sample_ID(i, :); % Extract sample ID
    brain_region = feather_to_mat_struct.Brain_Region(i, :); % Extract brain region
    time_series = transpose(feather_to_mat_struct.Time_Series(i, :)); % Extract time series

    % Compute fALFF for this participant/region
    SP_fALFF_res = SP_fALFF(time_series, samplingPeriod);
    
    fALFF = SP_fALFF_res.fALFF;
    ALFF = SP_fALFF_res.ALFF;

    % Add results to mat
    UCLA_CNP_fALFF_data{i, 1} = sample_ID;
    UCLA_CNP_fALFF_data{i, 2} = brain_region;
    UCLA_CNP_fALFF_data{i, 3} = fALFF;
    UCLA_CNP_fALFF_data{i, 4} = ALFF;
end

save(strcat(UCLA_CNP_data_path, 'processed_data/UCLA_CNP_AROMA_2P_GMR_ALFF_fALFF.mat'), 'UCLA_CNP_fALFF_data', '-v6')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ABIDE ASD
ABIDE_ASD_TS_file = strcat(ABIDE_ASD_data_path, "raw_data/ABIDE_ASD_FC1000_fMRI_TS.mat");
load(ABIDE_ASD_TS_file);

[m, n] = size(feather_to_mat_struct.Sample_ID);
% Initialise fALFF results
ABIDE_ASD_fALFF_data = cell(m,3);

for i = 1:m
    sample_ID = feather_to_mat_struct.Sample_ID(i, :); % Extract sample ID
    brain_region = feather_to_mat_struct.Brain_Region(i, :); % Extract brain region
    time_series = transpose(feather_to_mat_struct.Time_Series(i, :)); % Extract time series

    % Compute fALFF for this participant/region
    SP_fALFF_res = SP_fALFF(time_series, samplingPeriod);
    
    fALFF = SP_fALFF_res.fALFF;
    ALFF = SP_fALFF_res.ALFF;

    % Add results to mat
    ABIDE_ASD_fALFF_data{i, 1} = sample_ID;
    ABIDE_ASD_fALFF_data{i, 2} = brain_region;
    ABIDE_ASD_fALFF_data{i, 3} = fALFF;
    ABIDE_ASD_fALFF_data{i, 4} = ALFF;
end

save(strcat(ABIDE_ASD_data_path, 'processed_data/ABIDE_ASD_FC1000_ALFF_fALFF.mat'), 'ABIDE_ASD_fALFF_data', '-v6')
