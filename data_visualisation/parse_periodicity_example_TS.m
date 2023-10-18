% Load in the mat data
load("/Users/abry4213/data/TS_feature_manuscript/INP_Empirical1000.mat")

% Find index corresponding to MP_gingerbread_L300_IC_0.54_3.7_x
low_periodicity_TS_index = find(strcmp(labels, 'MP_gingerbread_L300_IC_0.54_3.7_x.dat'));
high_periodicity_TS_index = find(strcmp(labels, 'FL_duffvdp_L1000_N10000_IC_0.2_-0.25_y.dat'));

% Extract time series
low_periodicity_TS = timeSeriesData{low_periodicity_TS_index};
high_periodicity_TS = timeSeriesData{high_periodicity_TS_index};

% Create a timepoint variable for each array
low_timepoints = (1:length(low_periodicity_TS))';
high_timepoints = (1:length(high_periodicity_TS))';

% Combine into cell and write to CSV
low_table = table([low_periodicity_TS], ...
                      repmat({'Low'}, [length(low_periodicity_TS), 1]), ...
                      [low_timepoints], ...
                      'VariableNames', {'Value', 'Periodicity_Type', 'Timepoint'});


high_table = table([high_periodicity_TS], ...
                      repmat({'High'}, [length(high_periodicity_TS), 1]), ...
                      [high_timepoints], ...
                      'VariableNames', {'Value', 'Periodicity_Type', 'Timepoint'});

tables = vertcat(low_table, high_table);

% Write to CSV
writetable(tables, "/Users/abry4213/data/TS_feature_manuscript/periodicity_example_TS.csv")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Also compute the spline-detrended time series for two brain regions of
% interest
brain_region_TS = readtable("/Users/abry4213/data/TS_feature_manuscript/BOLD_data_for_wangs_PD.csv");
brain_region_TS_low = table2array(brain_region_TS(brain_region_TS.Brain_Region=="ctx-lh-caudalmiddlefrontal","values"));
brain_region_TS_high = table2array(brain_region_TS(brain_region_TS.Brain_Region=="ctx-lh-parstriangularis","values"));

% Create a timepoint variable for each array
low_timepoints = (1:length(brain_region_TS_low))';
high_timepoints = (1:length(brain_region_TS_high))';

% add catch22 files to path
addpath('/Users/abry4213/Library/CloudStorage/OneDrive-TheUniversityofSydney(Students)/github/catch22/Matlab/');

th = 0.01; % the threshold with which to count a peak

% Detrend the brain region with low catch22 periodicity
N = length(brain_region_TS_low); % length of the time series
ppSpline = splinefit(brain_region_TS_low'); % fit the three-knot cubic spline
y_spl = ppval(ppSpline,1:N); % find the fitted values
brain_region_TS_low_DT = brain_region_TS_low - y_spl'; % apply detrending

% Detrend the brain region with high catch22 periodicity
N = length(brain_region_TS_high); % length of the time series
ppSpline = splinefit(brain_region_TS_high'); % fit the three-knot cubic spline
y_spl = ppval(ppSpline,1:N); % find the fitted values
brain_region_TS_high_DT = brain_region_TS_high - y_spl'; % apply detrending

% Write the detrended time series to a table
low_table = table([brain_region_TS_low_DT], ...
                      repmat({'Low'}, [length(brain_region_TS_low_DT), 1]), ...
                      [low_timepoints], ...
                      'VariableNames', {'Detrended_Value', 'Periodicity_Type', 'Timepoint'});


high_table = table([brain_region_TS_high_DT], ...
                      repmat({'High'}, [length(brain_region_TS_high_DT), 1]), ...
                      [high_timepoints], ...
                      'VariableNames', {'Detrended_Value', 'Periodicity_Type', 'Timepoint'});

tables = vertcat(low_table, high_table);

% Write to CSV
writetable(tables, "/Users/abry4213/data/TS_feature_manuscript/brain_region_TS_detrended.csv")
