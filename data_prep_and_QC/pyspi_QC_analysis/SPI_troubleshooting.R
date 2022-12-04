################################################################################
# Load libraries
################################################################################

library(tidyverse)
library(icesTAF)
library(cowplot)
library(knitr)
library(kableExtra)
library(theft)
theme_set(theme_cowplot())

################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/QC/")
TAF::mkdir(plot_path)

SCZ_data_path <- "~/data/UCLA_Schizophrenia/"
SCZ_rdata_path <- paste0(SCZ_data_path, "processed_data/Rdata/")
SCZ_pydata_path <- paste0(SCZ_data_path, "raw_data/pydata/")
ASD_data_path <- "~/data/ABIDE_ASD/"
ASD_rdata_path <- paste0(ASD_data_path, "processed_data/Rdata/")
ASD_pydata_path <- paste0(ASD_data_path, "raw_data/pydata/")

# Define constants
SCZ_noise_proc <- "AROMA+2P+GMR"
ASD_noise_proc <- "FC1000"

# Load raw time-series data for UCLA Schizophrenia
SCZ_TS <- readRDS(paste0(SCZ_data_path, 
                         "raw_data/UCLA_Schizophrenia_fMRI_TS.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR")
# Load raw time-series data for ABIDE ASD
ASD_TS <- readRDS(paste0(ASD_data_path, 
                         "raw_data/ABIDE_ASD_fMRI_TS.Rds")) %>%
  filter(Noise_Proc == "FC1000")

# Load raw (non-normalized) pyspi14 data for UCLA Schizophrenia
SCZ_pyspi14 <- readRDS(paste0(SCZ_rdata_path, 
                              "UCLA_Schizophrenia_pyspi14_filtered.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR",
         Diagnosis %in% c("Control", "Schizophrenia"))
ASD_pyspi14 <- readRDS(paste0(ASD_rdata_path, 
                              "ABIDE_ASD_pyspi14_filtered.Rds")) %>%
  filter(Noise_Proc == "FC1000",
         Diagnosis %in% c("Control", "ASD"))

################################################################################
# Find SPIs that yielded NaN for at least one region pair and subject
# SCZ
SCZ_pyspi14 %>%
  filter(brain_region_from != brain_region_to) %>%
  filter(is.na(value)) %>%
  distinct(SPI) %>%
  pull(SPI)

res_round1 <- read.csv(paste0(SCZ_data_path, "raw_data/numpy_files/AROMA_2P_GMR/SPI_res_round1.csv")) %>%
  dplyr::select(-X) %>%
  filter(brain_region_from != brain_region_to) %>%
  dplyr::rename("value_1" = "value")
res_round2 <- read.csv(paste0(SCZ_data_path, "raw_data/numpy_files/AROMA_2P_GMR/SPI_res_round2.csv")) %>%
  dplyr::select(-X)%>%
  filter(brain_region_from != brain_region_to) %>%
  dplyr::rename("value_2" = "value")
merged = left_join(res_round1, res_round2)

diffs <- subset(merged, value_1 != value_2)

merged %>%
  filter(str_detect(SPI, "sgc")) %>%
  filter(is.nan(value_1) | is.nan(value_2))


to_follow_up_on <- SCZ_pyspi14 %>%
  filter(Sample_ID == "sub-10527",
         str_detect(SPI, "sgc"),
         brain_region_from != brain_region_to,
         is.nan(value)) %>%
  dplyr::select(brain_region_from, brain_region_to)
# ASD
ASD_pyspi14 %>%
  filter(brain_region_from != brain_region_to) %>%
  filter(is.na(value)) %>%
  distinct(SPI) %>%
  pull(SPI)

# Two SPIs yielded NaN for at least one region pair and subject:
# sgc_nonparametric_mean_fs-1_fmin-0_fmax-0-5
# di_gaussian

# Heatmap for each of these showing the subjects (rows)  and brain region pairs (columns)
plot_SPI_na <- function(this_SPI, 
                        SPI_name, 
                        pyspi_data) {
  p <- pyspi_data %>%
    filter(brain_region_from != brain_region_to,
           SPI == this_SPI) %>%
    mutate(fill = ifelse(is.na(value), T, F),
           region_pair = paste0(brain_region_from, brain_region_to)) %>%
    distinct(Sample_ID, region_pair, fill) %>%
    ggplot(data = ., mapping=aes(x=region_pair,
                                 y=Sample_ID,
                                 fill=fill)) +
    geom_tile() +
    scale_fill_manual(values = c("white", "red"), na.value="white") +
    ggtitle(paste0(this_SPI, "\nDistribution of NaN values")) +
    labs(fill = "NaN") +
    ylab("Subjects") +
    xlab("Region pairs") +
    theme(legend.position = "none",
          plot.title = element_text(hjust=0.5),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  return(p)
}
# UCLA Schizophrenia
plot_SPI_na(this_SPI="sgc_nonparametric_mean_fs-1_fmin-0_fmax-0-5",
            SPI_name="sgc_full_range",
            pyspi_data = SCZ_pyspi14)
ggsave(paste0(plot_path, "SCZ_sgc_full_range_NaN_heatmap_UCLA_Schizophrenia.png"),
       width=12, height=3, units="in", dpi=1200, bg="white")

plot_SPI_na(this_SPI="di_gaussian",
            SPI_name="di_gaussian",
            pyspi_data = SCZ_pyspi14)
ggsave(paste0(plot_path, SPI_name, "_NaN_heatmap_UCLA_Schizophrenia.png"),
       width=12, height=3, units="in", dpi=1200, bg="white")

# ABIDE ASD
plot_SPI_na(this_SPI="sgc_nonparametric_mean_fs-1_fmin-0_fmax-0-5",
            SPI_name="sgc_full_range",
            pyspi_data = ASD_pyspi14)
ggsave(paste0(plot_path, "sgc_full_range_NaN_heatmap_ABIDE_ASD.png"),
       width=12, height=8, units="in", dpi=300, bg="white")

plot_SPI_na(this_SPI="di_gaussian",
            SPI_name="di_gaussian",
            pyspi_data = ASD_pyspi14)
ggsave(paste0(plot_path, "di_gaussian_NaN_heatmap_ABIDE_ASD.png"),
       width=12, height=8, units="in", dpi=300, bg="white")



# First pick a region pair that yielded NaN for sgc full range from pyspi
# Focusing on sub-10527 with AROMA+2P+GMR noise processing
SCZ_TS %>%
  filter(Sample_ID == "sub-10527",
         Brain_Region %in% c("ctx-lh-caudalanteriorcingulate",
                             "ctx-lh-bankssts")) %>%
  mutate(Brain_Region = factor(Brain_Region,
                               levels = c("ctx-lh-caudalanteriorcingulate",
                                          "ctx-lh-bankssts"))) %>%
  ggplot(data=., mapping=aes(x=timepoint, y=values, color=Brain_Region)) +
  geom_line() +
  ggtitle("sgc_nonparametric_mean_fs-1_fmin-0_fmax-0-5\nthat yielded -0.498") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust=0.5)) +
  labs(color="Brain Region") +
  xlab("Time point") +
  ylab("AROMA+2P+GMR Signal Value") + 
  guides(color=guide_legend(nrow=2,byrow=TRUE))
ggsave(paste0(plot_path, "sgc_with_real_value_sub-10527_AROMA_2P_GMR.png"),
       width=6, height=4, units="in", dpi=300, bg="white")

# Plot left caudal anterior cingulate and left banks of the superior temporal
# sulcus for a subject for whom di_gaussian returned a real value (sub-10159) 
# versus a subject for whom di_gaussian returned NaN (sub-10171)
SCZ_TS %>%
  filter(Sample_ID %in% c("sub-10159", "sub-10171", "sub-50014"),
         Brain_Region %in% c("ctx-lh-caudalanteriorcingulate",
                             "ctx-lh-bankssts")) %>%
  rowwise() %>%
  mutate(Sample_Label = ifelse(Sample_ID == "sub-10159",
                               paste0(Sample_ID, ",\nReal"),
                               paste0(Sample_ID, ",\nNaN"))) %>%
  ungroup() %>%
  mutate(Brain_Region = factor(Brain_Region,
                               levels = c("ctx-lh-caudalanteriorcingulate",
                                          "ctx-lh-bankssts"))) %>%
  ggplot(data=., mapping=aes(x=timepoint, y=values, color=Brain_Region)) +
  geom_line() +
  facet_grid(Sample_Label ~ ., scales="fixed", switch="both") +
  ggtitle("Comparing raw fMRI time-series for di_gaussian") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust=0.5),
        strip.text.y.left = element_text(angle=0)) +
  labs(color="Brain Region") +
  xlab("Time point") +
  ylab("AROMA+2P+GMR Signal Value") + 
  guides(color=guide_legend(nrow=1,byrow=TRUE))
ggsave(paste0(plot_path, "di_gaussian_comparison_fMRI_TS_UCLA_Schizophrenia.png"),
       width=7, height=5, units="in", dpi=300, bg="white")


################################################################################

# I next tried re-running pyspi on the physics cluster with
# sgc_nonparametric_mean_fs-1_fmin-0-25_fmax-0-5 instead of sgc_nonparametric_mean_fs-1_fmin-0_fmax-0-5
# And I put di_gaussian at the top of the SPI computation list to see if
# Changing the order around made the difference

# I named this dataset with pyspi14_mod instead of pyspi14
SCZ_pyspi14_mod <- readRDS(paste0(SCZ_rdata_path, 
                               "UCLA_Schizophrenia_pyspi14_mod_filtered.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR")
ASD_pyspi14_mod <- readRDS(paste0(ASD_rdata_path, 
                                  "ABIDE_ASD_pyspi14_mod_filtered.Rds")) %>%
  filter(Noise_Proc == "FC1000")

################################################################################
# SPI robustness analysis
subset_to_subject_region <- function(input_TS_data,
                                     pydata_path,
                                     brain_regions,
                                     noise_proc,
                                     sample_ID,
                                     file_name) {
  input_TS_data %>%
    filter(Sample_ID == sample_ID,
           Noise_Proc == noise_proc,
           Brain_Region %in% brain_regions) %>%
    dplyr::select(Brain_Region, timepoint, values) %>%
    mutate(Brain_Region = factor(Brain_Region, levels = brain_regions)) %>%
    pivot_wider(id_cols = Brain_Region,
                names_from = "timepoint",
                values_from = "values") %>%
    arrange(Brain_Region) %>%
    dplyr::select(-Brain_Region) %>%
    write.table(., 
                file=paste0(pydata_path, file_name), 
                sep=",", col.names=F, row.names=F)
}

# Subject: sub-10159
# brain region 1: ctx-lh-bankssts
# brain region 2: ctx-lh-entorhinal
# noise processing: AROMA+2P+GMR
subset_to_subject_region(input_TS_data = SCZ_TS, pydata_path = SCZ_pydata_path,
                         brain_regions = c("ctx-lh-bankssts", "ctx-lh-entorhinal"),
                         noise_proc = SCZ_noise_proc,
                         sample_ID = "sub-10159",
                         file_name = "sub-10159_lh_bankssts_lh_entorhinal.csv")

# Subject: sub-10527
# brain region 1: ctx-lh-rostralanteriorcingulate
# brain region 2: ctx-lh-caudalmiddlefrontal
# noise processing: AROMA+2P+GMR
subset_to_subject_region(input_TS_data = SCZ_TS, pydata_path = SCZ_pydata_path,
                         brain_regions = c("ctx-lh-rostralanteriorcingulate", 
                                           "ctx-lh-caudalmiddlefrontal"),
                         noise_proc = SCZ_noise_proc,
                         sample_ID = "sub-10527",
                         file_name = "sub-10527_lh_rostralanteriorcingulate_lh_caudalmiddlefrontal.csv")

# Subject: 10021451277603445196
# brain region from: Precentral Gyrus
# brain region to: Angular Gyrus
# noise processing: AROMA+2P+GMR
subset_to_subject_region(input_TS_data = ASD_TS, pydata_path = ASD_pydata_path,
                         brain_regions = c("Precentral Gyrus",
                                           "Angular Gyrus"),
                         noise_proc = ASD_noise_proc,
                         sample_ID = "10021451277603445196",
                         file_name = "subject_10021451277603445196_precentral_angular.csv")

################################################################################
# Once data has been processed with pyspi 100x, visualise it here

all_SPI_100x_sub10159 <- read.csv(paste0(SCZ_pydata_path,
                                         "sub-10159_lh_bankssts_lh_entorhinal_all_SPIs.csv"),
                                  header=T) %>%
  dplyr::select(SPI, brain_region_from, brain_region_to, value, Iteration)

all_SPI_100x_sub10527 <- read.csv(paste0(SCZ_pydata_path,
                                         "sub-10527_lh_rostralanteriorcingulate_lh_caudalmiddlefrontal_all_SPIs.csv"),
                                  header=T) %>%
  dplyr::select(SPI, brain_region_from, brain_region_to, value, Iteration)

all_SPI_100x_sub10527_seeded <- read.csv(paste0(SCZ_pydata_path,
                                         "sub-10527_lh_rostralanteriorcingulate_lh_caudalmiddlefrontal_seeded_all_SPIs.csv"),
                                  header=T) %>%
  dplyr::select(SPI, brain_region_from, brain_region_to, value, Iteration)

all_SPI_100x_ASD_sub <- read.csv(paste0(ASD_pydata_path,
                                        "subject_10021451277603445196_precentral_angular_all_SPIs.csv"),
                                 header=T) %>%
  dplyr::select(SPI, brain_region_from, brain_region_to, value, Iteration)

# Function to find the non-deterministic SPIs in a dataset
find_nondeterministic_SPIs <- function(pyspi_res) {
  nd_SPIs <- pyspi_res %>%
    group_by(SPI, value) %>%
    count() %>%
    filter(n < 100) %>%
    ungroup() %>%
    group_by(SPI) %>%
    summarise(value_mean = mean(value, na.rm=T),
              value_SD = sd(value, na.rm=T),
              num_unique_values = n()) %>%
    arrange(desc(num_unique_values)) 
  
  return(nd_SPIs)
}

# Find SPIs where the value is not the same for all iterations
non_deterministic_SPIs_sub10159 <- find_nondeterministic_SPIs(all_SPI_100x_sub10159)
non_deterministic_SPIs_sub10527 <- find_nondeterministic_SPIs(all_SPI_100x_sub10527)
non_deterministic_SPIs_sub10527_seeded <- find_nondeterministic_SPIs(all_SPI_100x_sub10527_seeded)
non_deterministic_SPIs_ASD_sub <- find_nondeterministic_SPIs(all_SPI_100x_ASD_sub)

all_SPI_100x_ASD_sub %>%
  group_by(SPI, value) %>%
  mutate(n = n()) %>%
  filter(n < 100) %>%
  ungroup() %>%
  ggplot(data=., mapping=aes(x=value, fill=SPI)) +
  geom_histogram() +
  scale_fill_viridis_d()+ 
  xlab("SPI Value") +
  ylab("# Iterations") +
  facet_wrap(SPI ~ ., scales="free") +
  theme(legend.position = "none",
        axis.text.x = element_text(size=9))

# Find SPIs that yielded NaN
all_SPI_100x_sub10159 %>%
  filter(is.na(value)) %>%
  group_by(SPI) %>%
  count()%>%
  kable() %>%
  kable_styling(full_width = F)

all_SPI_100x_sub10527 %>%
  filter(is.na(value)) %>%
  group_by(SPI) %>%
  count()%>%
  kable() %>%
  kable_styling(full_width = F)

all_SPI_100x_ASD_sub %>%
  filter(is.na(value)) %>%
  group_by(SPI) %>%
  count()%>%
  kable() %>%
  kable_styling(full_width = F)