################################################################################
# Load libraries
################################################################################

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/opt/anaconda3/envs/pyspi/bin/python3"
reticulate::use_python(python_to_use)

library(tidyverse)
library(reticulate)
library(icesTAF)
library(cowplot)
library(ggpubr)
library(ggsignif)
library(patchwork)
library(feather)
library(glue)
library(R.matlab)
theme_set(theme_cowplot())

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
source(paste0(github_dir, "data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/movement_analysis/")
TAF::mkdir(plot_path)

study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Comparison_Group = c("Schizophrenia", "ADHD", "Bipolar", "ASD"),
                             Group_Title = c("SCZ", "ADHD", "BPD", "ASD"))

UCLA_CNP_data_path <- "~/data/UCLA_CNP/"
ABIDE_ASD_data_path <- "~/data/ABIDE_ASD/"

# Load data on subjects we actually used
UCLA_CNP_subjects_used <- pyarrow_feather$read_feather(glue("{UCLA_CNP_data_path}/processed_data/UCLA_CNP_filtered_sample_info_AROMA_2P_GMR_catch22_pyspi14.feather")) %>%
  pull(Sample_ID)
ABIDE_ASD_subjects_used <- pyarrow_feather$read_feather(glue("{ABIDE_ASD_data_path}/processed_data/ABIDE_ASD_filtered_sample_info_FC1000_catch22_pyspi14.feather")) %>%
  pull(Sample_ID)

# Load subject metadata
UCLA_CNP_sample_metadata <- feather::read_feather(glue("{UCLA_CNP_data_path}/study_metadata/UCLA_CNP_sample_metadata.feather")) %>%
  filter(Sample_ID %in% UCLA_CNP_subjects_used)
ABIDE_ASD_sample_metadata <- feather::read_feather(glue("{ABIDE_ASD_data_path}/study_metadata/ABIDE_ASD_sample_metadata.feather")) %>%
  filter(Sample_ID %in% ABIDE_ASD_subjects_used)

# Load mean framewise displacement data
UCLA_CNP_mean_FD <- read.table(glue("{UCLA_CNP_data_path}/movement_data/fmriprep/UCLA_CNP_mFD.txt"), 
                                       sep=",", colClasses = "character")
ABIDE_ASD_mean_FD <- read.table(glue("{ABIDE_ASD_data_path}/movement_data/fmriprep/ABIDE_ASD_mFD.txt"), 
                                     sep=",", colClasses = "character")
# Load full framewise displacement data
UCLA_CNP_full_FD <- as.data.frame(readMat(glue("{UCLA_CNP_data_path}/movement_data/fmriprep/UCLA_CNP_all_FD.mat"))[[1]])
ABIDE_ASD_full_FD <- as.data.frame(readMat(glue("{ABIDE_ASD_data_path}/movement_data/fmriprep/ABIDE_ASD_all_FD.mat"))[[1]])
colnames(UCLA_CNP_full_FD) <- colnames(ABIDE_ASD_full_FD) <- colnames(UCLA_CNP_mean_FD) <- colnames(ABIDE_ASD_mean_FD) <- c("Sample_ID", "Jenkinson", "Power", "VanDijk")

# Un-list full framewise displacement data
UCLA_CNP_full_FD <- UCLA_CNP_full_FD %>%
  unnest(cols=c(Sample_ID, Jenkinson, Power,VanDijk)) %>%
  unnest(cols=c(Sample_ID, Jenkinson, Power, VanDijk)) %>%
  filter(Sample_ID %in% UCLA_CNP_subjects_used) %>%
  group_by(Sample_ID) %>%
  mutate(Frame_Number = row_number())
ABIDE_ASD_full_FD <- ABIDE_ASD_full_FD %>%
  unnest(cols=c(Sample_ID, Jenkinson, Power,VanDijk)) %>%
  unnest(cols=c(Sample_ID, Jenkinson, Power, VanDijk)) %>%
  filter(Sample_ID %in% ABIDE_ASD_subjects_used) %>%
  group_by(Sample_ID) %>%
  mutate(Frame_Number = row_number())
colnames(UCLA_CNP_full_FD) <- colnames(ABIDE_ASD_full_FD) <- c("Sample_ID", "Jenkinson", "Power", "VanDijk", "Frame_Number")
  
# Set mFD columns as numeric
UCLA_CNP_mean_FD <- UCLA_CNP_mean_FD %>%
  mutate_at(c("Jenkinson", "Power", "VanDijk"), function(x) as.numeric(x)) %>%
  left_join(., UCLA_CNP_sample_metadata) %>%
  filter(Sample_ID %in% UCLA_CNP_subjects_used)
ABIDE_ASD_mean_FD <- ABIDE_ASD_mean_FD %>%
  mutate_at(c("Jenkinson", "Power", "VanDijk"), function(x) as.numeric(x)) %>%
  left_join(., ABIDE_ASD_sample_metadata) %>%
  filter(Sample_ID %in% ABIDE_ASD_subjects_used)

################################################################################
# Compare FD-Power distributions between each case-control comparison
################################################################################

plot_group_vs_control <- function(FD_dataset,
                                  study,
                                  dx,
                                  dx_title) {
  p <- FD_dataset %>%
    filter(Study==study,
           Diagnosis %in% c("Control", dx)) %>%
    mutate(Diagnosis = factor(Diagnosis, levels = c(dx, "Control"))) %>%
    ggplot(data=., mapping=aes(x=Diagnosis, y=Power)) +
    geom_violin(aes(fill=Diagnosis)) +
    geom_boxplot(color="black", fill=NA, width=0.1) +
    geom_signif(test = "wilcox.test",
                comparisons = list(c(dx, "Control")), 
                map_signif_level=TRUE) +
    scale_fill_manual(values = c("#00B06D", "#737373")) +
    scale_y_continuous(expand = c(0,0,0.1,0)) +
    ggtitle(dx_title) +
    ylab("mFD-Power") +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          plot.title = element_text(hjust=0.5))
  return(p)
}

plots <- 1:4 %>%
  purrr::map(~ plot_group_vs_control(plyr::rbind.fill(UCLA_CNP_mean_FD, ABIDE_ASD_mean_FD),
                                     study = study_group_df$Study[.x],
                                     dx = study_group_df$Comparison_Group[.x],
                                     dx_title = study_group_df$Group_Title[.x]))

wrap_plots(plots, ncol=2)
ggsave(paste0(plot_path, "mFD_Power_by_Group.png"),
       width = 6, height=5, units="in", dpi=300)
