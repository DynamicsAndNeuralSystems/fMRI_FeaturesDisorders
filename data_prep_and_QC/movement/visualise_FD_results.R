################################################################################
# Script to compare framewise displacement (FD) using Jenkinson et al. (2002) 
# root-mean-square (RMS) approach, comparing Linden's results with Annie's
# results

# Author: Annie G. Bryant, 28 November 2022
################################################################################

# Load libraries
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

################################################################################
# UCLA Schizophrenia dataset
################################################################################

# Define where data is stored
SCZ_movement_data_path = "~/data/UCLA_Schizophrenia/movementData/"
plot_path = "~/github/fMRI_FeaturesDisorders/plots/QC/"

# Read in Linden's mean FD results
SCZ_linden_FD_res <- read.table(paste0(SCZ_movement_data_path, "fdAvgs_UCLA_Schizophrenia.txt"))
colnames(SCZ_linden_FD_res) <- "mFD_Linden"

# Read in Annie's mean FD results
SCZ_annie_FD_res <- read.table(paste0(SCZ_movement_data_path, "fdAvgs_UCLA_Schizophrenia_Annie.txt"),
                           sep=",")
colnames(SCZ_annie_FD_res) <- c("Sample_ID", "mFD_Annie")

# Merge Linden's results with Annie's results
SCZ_merged_FD_res <- cbind(SCZ_annie_FD_res, SCZ_linden_FD_res)

# Plot the correlation between Linden and Annie mean FDjenk calculations
SCZ_merged_FD_res %>%
  ggplot(data=., mapping=aes(x=mFD_Linden,  y=mFD_Annie)) +
  geom_point() +
  ggtitle("Framewise Displacement (FD) with\nJenkinson 2002 RMS Method") +
  ylab("Mean FDjenk (Annie)") +
  xlab("Mean FDjenk (Linden)") +
  geom_abline(color="red") +
  coord_equal() +
  theme(plot.title = element_text(hjust=0.5))
ggsave(paste0(plot_path, "Movement_FD_Linden_vs_Annie_UCLA_Schizophrenia.png"),
       width=5, height=4, dpi=300, bg="white")

# Plot distribution of FDjenk values with Annie's results
SCZ_annie_FD_res %>%
  ggplot(data=., mapping=aes(x=mFD_Annie)) +
  geom_histogram(bins=50, fill="lightsteelblue") +
  ggtitle("FDjenk Distribution in\nUCLA Schizophrenia Dataset") +
  xlab("Mean FDjenk (Annie)") +
  ylab("# Participants") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(paste0(plot_path, "Movement_FDjenk_Annie_dist_UCLA_Schizophrenia.png"),
       width=4.5, height=4, dpi=300, bg="white")


################################################################################
# ABIDE ASD dataset
################################################################################

# Define where data is stored
ASD_movement_data_path = "~/data/ABIDE_ASD/movement_data/"
plot_path = "~/github/fMRI_FeaturesDisorders/plots/QC/"

# Read in Annie's mean FD results
ASD_annie_FD_res <- read.table(paste0(ASD_movement_data_path, "fdAvgs_ABIDE_ASD_Annie.txt"),
                               sep=",")
colnames(ASD_annie_FD_res) <- c("Sample_ID", "mFD_Annie")

# Plot distribution of FDjenk values with Annie's results
ASD_annie_FD_res %>%
  ggplot(data=., mapping=aes(x=mFD_Annie)) +
  geom_histogram(bins=50, fill="lightsteelblue") +
  ggtitle("FDjenk Distribution in\nABIDE ASD Dataset") +
  xlab("Mean FDjenk (Annie)") +
  ylab("# Participants") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(paste0(plot_path, "Movement_FDjenk_Annie_dist_ABIDE_ASD.png"),
       width=4.5, height=4, dpi=300, bg="white")

# Divide mean FDjenk by 100 and re-plot
ASD_annie_FD_res %>%
  ggplot(data=., mapping=aes(x=mFD_Annie/100)) +
  geom_histogram(bins=50, fill="lightsteelblue") +
  ggtitle("0.01*FDjenk Distribution in\nABIDE ASD Dataset") +
  xlab("Mean FDjenk (Annie) divided by 100") +
  ylab("# Participants") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(paste0(plot_path, "Movement_FDjenk_Annie_dist_ABIDE_ASD_Divided.png"),
       width=4.5, height=4, dpi=300, bg="white")
