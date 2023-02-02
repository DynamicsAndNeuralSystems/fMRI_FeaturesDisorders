% Calculate mean framewise displacement (FD) using root-mean-square
% approach outlined in Jenkinson et al. (NeuroImage) 2002

% Functions used from github.com/lindenparkes/rsfMRI/


%%%%%%%%%%%%%%%%%%%%%%%%%%%% UCLA CNP %%%%%%%%%%%%%%%%%%%%%%%%%%%
UCLA_CNP_subject_movement_data_path = "/headnode1/abry4213/data/UCLA_CNP/movement_data/";
txt_files = dir(fullfile(UCLA_CNP_subject_movement_data_path, "*movData.txt"));

% Instantiate matrix with two columns, one for subject ID and one for mean
% FD
num_subjects = length(txt_files);
FD_m_mat = repmat("hello", num_subjects, 4);

for i = 1:num_subjects
    % Define subject ID and movement data file
    file = txt_files(i).name;
    subject = regexprep(file, "_movData.txt", "");

    % Read in movement data
    subject_movement_data = dlmread(fullfile(UCLA_CNP_subject_movement_data_path, file));

    % Calculate framewise displacement
    % Jenkinson formula
    subject_FD_jenk = GetFDJenk(subject_movement_data);
    % Power formula
    subject_FD_power = GetFDPower(subject_movement_data);
    % VanD formula
    subject_FD_VanD = GetFDVanD(subject_movement_data);

    % Calculate and store mean FD
    subject_m_FD_jenk = mean(subject_FD_jenk);
    subject_m_FD_power = mean(subject_FD_power);
    subject_m_FD_VanD = mean(subject_FD_VanD);
    FD_m_mat(i,1) = subject;
    FD_m_mat(i,2) = subject_m_FD_jenk;
    FD_m_mat(i,3) = subject_m_FD_power;
    FD_m_mat(i,4) = subject_m_FD_VanD;
end

% Write my results to a .txt file
output_file = fullfile(UCLA_CNP_subject_movement_data_path, "UCLA_CNP_mFD.txt");
writematrix(FD_m_mat, output_file)
