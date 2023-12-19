function exitcode = compute_power_FD(dataset_ID, github_path, movement_data_path)
    % Calculate mean framewise displacement (FD) using root-mean-square
    % approach outlined in Jenkinson et al. (NeuroImage) 2002
    
    % Functions used from github.com/lindenparkes/rsfMRI/
    
    % Read in movement data
    txt_files = dir(fullfile(movement_data_path, "*movData.txt"));
    
    % Add data path to path
    addpath(genpath(github_path));
    % addpath(movement_data_path);
    cd(movement_data_path);
    
    % Instantiate matrix with two columns, one for subject ID and one for mean
    % FD
    num_subjects = length(txt_files);
    FD_mat = num2cell(repmat("hello", num_subjects, 2));
    FD_m_mat = repmat("hello", num_subjects, 2);
    
    for i = 1:num_subjects
        % Define subject ID and movement data file
        file = txt_files(i).name;
        subject = regexprep(file, "_movData.txt", "");
    
        % Read in movement data
        subject_movement_data = dlmread(file);
    
        % NOTE: ABIDE ASD data was processed with AFNI rather than SPM
        % 3dvolreg outputs (rot,trans) rather than SPM (trans,rot)
        if dataset_ID == "ABIDE"
            subject_movement_data = subject_movement_data(:, [4:6, 1:3]);
        end
    
        % Calculate framewise displacement with Power formula
        subject_FD_power = GetFDPower(subject_movement_data);
    
        % Calculate and store mean FD
        subject_m_FD_power = mean(subject_FD_power);
        FD_m_mat(i,1) = subject;
        FD_m_mat(i,2) = subject_m_FD_power;
    end
    
    % Write mean FD results to a .txt file for easy viewing
    output_mean_file = strcat(dataset_ID, '_Mean_FD_Power.txt');
    writematrix(FD_m_mat, output_mean_file);

end