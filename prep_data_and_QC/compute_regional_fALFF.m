function []=prog(data_path, SP_fALFF_script_path, TS_mat_file, output_mat_file)

% Add SP_fALFF script to path
addpath(SP_fALFF_script_path);

% Define TR as 2s
samplingPeriod = 2;

% Load in feather data file
load(TS_mat_file);

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
end

save(output_mat_file, 'fALFF_data', '-v6')
