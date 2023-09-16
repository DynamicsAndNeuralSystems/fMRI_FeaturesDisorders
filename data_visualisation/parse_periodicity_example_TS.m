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