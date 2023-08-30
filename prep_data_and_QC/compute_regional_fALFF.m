UCLA_CNP_data_path = '/Users/abry4213/data/UCLA_CNP/raw_data/';
ABIDE_ASD_data_path = '/Users/abry4213/data/ABIDE_ASD/raw_data/';

% Load in feather data file
UCLA_CNP_TS_file = "UCLA_CNP_AROMA_2P_GMR_fMRI_TS.feather";

%%
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

%% 

ALFF_data = cell(size(TimeSeries,1),3);
samplingPeriod = 2.68;

for i = 1:size(TimeSeries,1)
    TS_sample = TimeSeries.Name{i};
    TS_data = TimeSeries.Data{i};

    % Ben's implementation of fALFF in Fallon et al 2020
    fALFF = SP_fALFF(TS_data, samplingPeriod);
    % Append to ALFF_data
    ALFF_data{i, 1} = TS_sample;
    ALFF_data{i, 2} = fALFF.ALFF;
    ALFF_data{i, 3} = fALFF.fALFF;
end

save(strcat(data_path, 'All_rsfMRI_data_ALFF_fALFF.mat'), 'ALFF_data', '-v6')