

#-------------------------------------------------------------------------------
# Function to read in file with average fractional displacement (FD)
# As well as list of individual subject movement files
# And output a CSV containing the Subject ID, diagnosis, and FD
#-------------------------------------------------------------------------------
compile_movement_data <- function(fd_path, 
                                  input_dataset_name,
                                  subject_csv) {
  mov_data <- read.table(paste0(fd_path, sprintf("fdAvgs_%s.txt",
                                                 input_dataset_name)))
  colnames(mov_data)[1] <- "FD"
  
  mov_data_subjects <- list.files(fd_path, pattern="_movData.txt") %>%
    gsub("_movData.txt", "", .)
  
  mov_data$Subject_ID <- mov_data_subjects
  
  subject_info <- read.csv(subject_csv)
  colnames(subject_info)[1] <- "Subject_ID"
  
  mov_data <- mov_data %>%
    left_join(., subject_info) %>%
    dplyr::select(Subject_ID, diagnosis, FD)
  
  return(mov_data)
}

#-------------------------------------------------------------------------------
# Plot FD values by diagnosis boxplot
#-------------------------------------------------------------------------------
plot_FD_vs_diagnosis <- function(movement_data, input_dataset_name) {
  movement_data %>%
    ggplot(data=., mapping=aes(x=diagnosis, y=FD, fill=diagnosis)) +
    geom_boxplot() +
    ylab("Fractional Displacement (FD)") +
    xlab("Diagnosis") +
    ggtitle(sprintf("%s Subject Movement by Diagnosis",
                    input_dataset_name)) +
    theme(legend.position="none",
          plot.title=element_text(hjust=0.5)) 
}

#-------------------------------------------------------------------------------
# Plot the number of control vs. schizophrenia subjects retained per FD threshold
#-------------------------------------------------------------------------------
plot_subjects_per_fd_threshold <- function(movement_data,
                                           input_dataset_name) {
  fd_thresh_list <- list()
  for (fd_threshold in seq(0, 1, by=0.01)) {
    num_ctrl = nrow(subset(movement_data, FD <= fd_threshold & diagnosis=="CONTROL"))
    num_schz = nrow(subset(movement_data, FD <= fd_threshold & diagnosis=="SCHZ"))
    thresh_df <- data.frame(FD_Threshold = fd_threshold,
                            Control = num_ctrl,
                            Schizophrenia = num_schz)
    fd_thresh_list <- rlist::list.append(fd_thresh_list, thresh_df)
  }
  threshold_data <- do.call(plyr::rbind.fill, fd_thresh_list)
  
  threshold_data %>%
    pivot_longer(cols=c(-FD_Threshold),
                 names_to="Group",
                 values_to="# Subjects") %>%
    group_by(Group) %>%
    mutate(`% Subjects` = `# Subjects` / max(`# Subjects`, na.rm=T)) %>%
    pivot_longer(cols = c(`# Subjects`, `% Subjects`)) %>%
    ggplot(data=., mapping=aes(x=FD_Threshold, y=value, color=Group, group=Group)) +
    geom_line(size=2) +
    xlab("FD Threshold") +
    scale_x_reverse() +
    facet_grid(name ~ ., scales="free", switch="both") +
    ggtitle(sprintf("%s Subjects\nRetained by FD Threshold",
                    input_dataset_name)) +
    theme(legend.position="bottom",
          strip.placement = "outside",
          plot.title=element_text(hjust=0.5))
}

#-------------------------------------------------------------------------------
# Plot the ratio of schizophrenia:control subjects retained per FD threshold
#-------------------------------------------------------------------------------

plot_schz_ctrl_ratio_per_fd_threshold <- function(movement_data,
                                                  input_dataset_name) {
  fd_thresh_list <- list()
  for (fd_threshold in seq(0, 1, by=0.01)) {
    num_ctrl = nrow(subset(movement_data, 
                           FD <= fd_threshold & diagnosis=="CONTROL"))
    num_schz = nrow(subset(movement_data, 
                           FD <= fd_threshold & diagnosis=="SCHZ"))
    thresh_df <- data.frame(FD_Threshold = fd_threshold,
                            Control = num_ctrl,
                            Schizophrenia = num_schz)
    fd_thresh_list <- rlist::list.append(fd_thresh_list, thresh_df)
  }
  threshold_data <- do.call(plyr::rbind.fill, fd_thresh_list) %>%
    mutate(SCZ_to_CTRL = Schizophrenia / Control)
  
  threshold_data %>%
    ggplot(data=., mapping=aes(x=FD_Threshold, y=SCZ_to_CTRL)) +
    geom_line(size=2) +
    ylab("SCZ:CTRL Ratio") +
    xlab("FD Threshold") +
    scale_x_reverse() +
    ggtitle(sprintf("Ratio of SCZ:CTRL %s Subjects\nRetained by FD Threshold",
                    input_dataset_name)) +
    theme(legend.position="bottom",
          plot.title=element_text(hjust=0.5))
}