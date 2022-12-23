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
colnames(SCZ_annie_FD_res) <- c("Sample_ID", "mFD_Annie_Jenk", "mFD_Annie_Power", "mFD_Annie_VanDijk")

# Merge Linden's results with Annie's results
SCZ_merged_FD_res <- cbind(SCZ_annie_FD_res, SCZ_linden_FD_res)

# Plot the correlation between Linden and Annie mean FDjenk calculations
SCZ_merged_FD_res %>%
  ggplot(data=., mapping=aes(x=mFD_Linden,  y=mFD_Annie_Jenk)) +
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
  pivot_longer(cols=c(-Sample_ID),
               names_to="Metric",
               values_to = "Value") %>%
  mutate(Metric = gsub("mFD_Annie_", "", Metric)) %>%
  ggplot(data=., mapping=aes(x=Value)) +
  geom_histogram(bins=50, fill="lightsteelblue") +
  ggtitle("FD Distribution in\nUCLA Schizophrenia Dataset") +
  xlab("Mean FD (Annie)") +
  ylab("# Participants") +
  facet_wrap(Metric ~ ., scales="free") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(paste0(plot_path, "Movement_FD_Annie_dist_UCLA_Schizophrenia.png"),
       width=7, height=3, dpi=300, bg="white")


################################################################################
# ABIDE ASD dataset
################################################################################

# Define where data is stored
ASD_movement_data_path = "~/data/ABIDE_ASD/movement_data/"
plot_path = "~/github/fMRI_FeaturesDisorders/plots/QC/"

# Read in Annie's mean FD results
ASD_annie_FD_res <- read.table(paste0(ASD_movement_data_path, "fdAvgs_ABIDE_ASD_Annie.txt"),
                               sep=",")
colnames(ASD_annie_FD_res) <- c("Sample_ID", "mFD_Annie_Jenk", "mFD_Annie_Power", "mFD_Annie_VanDijk")

ASD_annie_FD_res %>%
  pivot_longer(cols=c(-Sample_ID),
               names_to="Metric",
               values_to = "Value") %>%
  mutate(Metric = gsub("mFD_Annie_", "", Metric)) %>%
  ggplot(data=., mapping=aes(x=Value)) +
  geom_histogram(bins=50, fill="lightsteelblue") +
  ggtitle("FD Distribution in\nABIDE ASD Dataset") +
  xlab("Mean FD (Annie)") +
  ylab("# Participants") +
  facet_wrap(Metric ~ ., scales="free") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(paste0(plot_path, "Movement_FD_Annie_dist_ABIDE_ASD.png"),
       width=7, height=3, dpi=300, bg="white")

# Filter to FD < 5
ASD_annie_FD_res %>%
  pivot_longer(cols=c(-Sample_ID),
               names_to="Metric",
               values_to = "Value") %>%
  filter(Value < 5) %>%
  mutate(Metric = gsub("mFD_Annie_", "", Metric)) %>%
  ggplot(data=., mapping=aes(x=Value)) +
  geom_histogram(bins=50, fill="lightsteelblue") +
  ggtitle("FD Distribution in\nABIDE ASD Dataset") +
  xlab("Mean FD (Annie)") +
  ylab("# Participants") +
  facet_wrap(Metric ~ ., scales="free") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(paste0(plot_path, "Movement_FD_Annie_dist_ABIDE_ASD_under5.png"),
       width=7, height=3, dpi=300, bg="white")


summary(ASD_annie_FD_res)
