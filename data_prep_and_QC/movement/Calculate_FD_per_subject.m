% Calculate mean framewise displacement (FD) using root-mean-square
% approach outlined in Jenkinson et al. (NeuroImage) 2002

% Functions used from github.com/lindenparkes/rsfMRI/


%%%%%%%%%%%%%%%%%%%%%%%%%%%% UCLA Schizophrenia %%%%%%%%%%%%%%%%%%%%%%%%%%%
SCZ_subject_movement_data_path = "/Users/abry4213/data/UCLA_Schizophrenia/movementData/";
txt_files = dir(fullfile(SCZ_subject_movement_data_path, "*movData.txt"));

% Instantiate matrix with two columns, one for subject ID and one for mean
% FD
num_subjects = length(txt_files);
FD_m_mat = repmat("hello", num_subjects, 2);

for i = 1:num_subjects
    % Define subject ID and movement data file
    file = txt_files(i).name;
    subject = regexprep(file, "_movData.txt", "");

    % Read in movement data
    subject_movement_data = dlmread(fullfile(SCZ_subject_movement_data_path, file));

    % Calculate framewise displacement with Jenkinson formula
    subject_FD = GetFDJenk(subject_movement_data);

    % Calculate and store mean FD
    subject_m_FD = mean(subject_FD);
    FD_m_mat(i,1) = subject;
    FD_m_mat(i,2) = subject_m_FD;
end

% Write my results to a .txt file
output_file = fullfile(SCZ_subject_movement_data_path, "fdAvgs_UCLA_Schizophrenia_Annie.txt");
writematrix(FD_m_mat, output_file)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ABIDE ASD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ASD_subject_movement_data_path = "/Users/abry4213/data/ABIDE_ASD/movement_data/";
txt_files = dir(fullfile(ASD_subject_movement_data_path, "*movData.txt"));

% Instantiate matrix with two columns, one for subject ID and one for mean
% FD
num_subjects = length(txt_files);
FD_m_mat = repmat("hello", num_subjects, 2);

for i = 1:num_subjects
    % Define subject ID and movement data file
    file = txt_files(i).name;
    subject = regexprep(file, "_movData.txt", "");

    % Read in movement data
    subject_movement_data = dlmread(fullfile(ASD_subject_movement_data_path, file));
    % NOTE: ABDIE ASD data was processed with AFNI rather than SPM
    % 3dvolreg outputs (rot,trans) rather than SPM's (trans,rot)
    subject_movement_data_reordered = subject_movement_data(:, [4:6, 1:3]);

    % Convert degrees to radians for rotation columns
    subject_movement_data_reordered_degrees = subject_movement_data_reordered;
    subject_movement_data_reordered_degrees(:, [4:6]) = deg2rad(subject_movement_data_reordered(:, [4:6]));

    % Calculate framewise displacement with Jenkinson formula
    subject_FD = GetFDJenk(subject_movement_data_reordered);

    % Calculate and store mean FD
    subject_m_FD = mean(subject_FD);
    FD_m_mat(i,1) = subject;
    FD_m_mat(i,2) = subject_m_FD;
end

% Write my results to a .txt file
output_file = fullfile(ASD_subject_movement_data_path, "fdAvgs_ABIDE_ASD_Annie.txt");
writematrix(FD_m_mat, output_file)