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
ASD_data_path <- "~/data/ABIDE_ASD/"
ASD_rdata_path <- paste0(ASD_data_path, "processed_data/Rdata/")

# Load raw (non-normalized) pyspi14 data for UCLA Schizophrenia
SCZ_pyspi14 <- readRDS(paste0(SCZ_rdata_path, 
                              "UCLA_Schizophrenia_pyspi14_filtered.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR",
         Diagnosis %in% c("Control", "Schizophrenia"))
ASD_pyspi14 <- readRDS(paste0(ASD_rdata_path, 
                              "ABIDE_ASD_pyspi14_filtered.Rds")) %>%
  filter(Noise_Proc == "FC1000",
         Diagnosis %in% c("Control", "ASD"))

# Find SPIs that yielded NaN for at least one region pair and subject
# SCZ
SCZ_pyspi14 %>%
  filter(brain_region_1 != brain_region_2) %>%
  filter(is.na(value)) %>%
  distinct(SPI) %>%
  pull(SPI)
# ASD
ASD_pyspi14 %>%
  filter(brain_region_1 != brain_region_2) %>%
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
    filter(brain_region_1 != brain_region_2,
           SPI == this_SPI) %>%
    mutate(fill = ifelse(is.na(value), T, F),
           region_pair = paste0(brain_region_1, brain_region_2)) %>%
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

# Load raw time-series data for UCLA Schizophrenia
SCZ_TS <- readRDS(paste0(SCZ_data_path, 
                          "raw_data/UCLA_Schizophrenia_fMRI_TS.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR")

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

# Check distribution of remaining NA for SGC
### UCLA Schizophrenia
SCZ_pyspi14_mod %>%
  filter(SPI=="sgc_nonparametric_mean_fs-1_fmin-0-25_fmax-0-5",
         brain_region_1 != brain_region_2) %>%
  filter(is.na(value)) %>%
  arrange(brain_region_1, brain_region_2, Sample_ID) %>%
  dplyr::select(-SPI, -Noise_Proc) %>%
  kable(.) %>%
  kable_styling(full_width=F)
# Original full-range SGC comparison
SCZ_pyspi14 %>%
  filter(SPI=="sgc_nonparametric_mean_fs-1_fmin-0_fmax-0-5",
         brain_region_1 != brain_region_2) %>%
  filter(is.na(value)) %>%
  arrange(brain_region_1, brain_region_2, Sample_ID) %>%
  dplyr::select(-SPI, -Noise_Proc) 

# Find # NaN for SGC for each brain region
SCZ_pyspi14_mod %>%
  filter(SPI=="sgc_nonparametric_mean_fs-1_fmin-0-25_fmax-0-5",
         brain_region_1 != brain_region_2) %>%
  filter(is.na(value)) %>%
  dplyr::select(Sample_ID, brain_region_1, brain_region_2) %>%
  pivot_longer(cols=c(-Sample_ID),
               values_to="brain_region") %>%
  mutate(hemisphere = ifelse(str_detect(brain_region, "rh-|Right-"),
                             "Right", "Left")) %>%
  mutate(brain_region = gsub("ctx-rh-|ctx-lh-|Right-|Left-", "", brain_region)) %>%
  group_by(brain_region, hemisphere) %>%
  count() %>%
  mutate(n = ifelse(hemisphere=="Left", -1*n, n)) %>%
  mutate(hemisphere = factor(hemisphere, levels = c("Left", "Right"))) %>%
  ggplot(data=., mapping=aes(x=n, y=fct_reorder(brain_region,
                                                abs(n),
                                                .fun=max,
                                                .desc=F), fill=hemisphere)) +
  geom_col() +
  scale_x_continuous(breaks = seq(-4, 4, by=2),
                     labels = abs(seq(-4, 4, by=2))) +
  ggtitle("sgc_nonparametric_mean_fs-1_fmin-0-25_fmax-0-5\nRemaining NaNs by Brain Region") +
  labs(fill = "Hemisphere") +
  ylab("Brain Region") +
  xlab("# NaN in/out") +
  theme(legend.position="bottom",
        plot.title=element_text(hjust=0.5))
ggsave(paste0(plot_path, "sgc_nonparametric_mean_fs-1_fmin-0-25_fmax-0-5_NaN_UCLA_Schizophrenia_pyramid.png"),
       width=8, height=5, units="in", dpi=300, bg="white")

### ABIDE ASD
ASD_pyspi14_mod %>%
  filter(SPI=="sgc_nonparametric_mean_fs-1_fmin-0-25_fmax-0-5",
         brain_region_1 != brain_region_2) %>%
  filter(is.na(value)) %>%
  arrange(brain_region_1, brain_region_2, Sample_ID) %>%
  dplyr::select(-SPI, -Noise_Proc)
# Original full-range SGC comparison
ASD_pyspi14 %>%
  filter(SPI=="sgc_nonparametric_mean_fs-1_fmin-0_fmax-0-5",
         brain_region_1 != brain_region_2) %>%
  filter(is.na(value)) %>%
  arrange(brain_region_1, brain_region_2, Sample_ID) %>%
  dplyr::select(-SPI, -Noise_Proc)  %>%
  distinct(Sample_ID)
# Find subjects with NaN still for 0.25-0.5 range
ASD_pyspi14_mod %>%
  filter(SPI=="sgc_nonparametric_mean_fs-1_fmin-0-25_fmax-0-5",
         brain_region_1 != brain_region_2) %>%
  filter(is.na(value)) %>%
  arrange(brain_region_1, brain_region_2, Sample_ID) %>%
  dplyr::select(-SPI, -Noise_Proc)

# Find # NaN for SGC for each brain region
ASD_pyspi14_mod %>%
  filter(SPI=="sgc_nonparametric_mean_fs-1_fmin-0-25_fmax-0-5",
         brain_region_1 != brain_region_2) %>%
  filter(is.na(value)) %>%
  dplyr::select(Sample_ID, brain_region_1, brain_region_2) %>%
  pivot_longer(cols=c(-Sample_ID),
               values_to="brain_region") %>%
  group_by(brain_region) %>%
  count() %>%
  ggplot(data=., mapping=aes(y=n, x=fct_reorder(brain_region,
                                                n,
                                                .fun=max,
                                                .desc=T))) +
  geom_col(aes(fill = n)) +
  ggtitle("sgc_nonparametric_mean_fs-1_fmin-0-25_fmax-0-5\nRemaining NaNs by Brain Region") +
  xlab("Brain Region") +
  ylab("# NaN in/out") + 
  scale_fill_viridis_c() +
  scale_x_discrete(labels = scales::wrap_format(30))+
  theme(legend.position="none",
        axis.text.x = element_text(angle=90, hjust=1, size=8, vjust=0.4),
        plot.title=element_text(hjust=0.5))
ggsave(paste0(plot_path, "sgc_nonparametric_mean_fs-1_fmin-0-25_fmax-0-5_NaN_ABIDE_ASD.png"),
       width=15, height=5, units="in", dpi=300, bg="white")

# di_gaussian redo
  
### UCLA Schizophrenia
SCZ_pyspi14_mod %>%
  filter(SPI=="di_gaussian",
         brain_region_1 != brain_region_2) %>%
  filter(is.na(value)) %>%
  group_by(Sample_ID, Diagnosis) %>%
  summarise(num_NA = n()) %>%
  ungroup() %>%
  arrange(desc(num_NA)) %>%
  kable(.) %>%
  kable_styling(full_width=F)

plot_SPI_na(this_SPI="di_gaussian",
            SPI_name="di_gaussian_mod",
            pyspi_data = SCZ_pyspi14_mod)

SCZ_pyspi14_mod %>%
  filter(SPI=="di_gaussian",
         brain_region_1 != brain_region_2) %>%
  filter(is.na(value)) %>%
  group_by(Sample_ID, Diagnosis) %>%
  summarise(num_NA = n()) %>%
  ungroup() %>%
  pull(Sample_ID)

# I re-ran just di_gaussian for these 13 subjects with one or more NaN directly 
# within spyder on the physics cluster
SCZ_di_gauss_mod_NaN <- SCZ_pyspi14_mod %>%
  filter(SPI=="di_gaussian",
         brain_region_1 != brain_region_2) %>%
  filter(is.na(value)) %>%
  group_by(Sample_ID, Diagnosis) %>%
  summarise(num_NA = n()) %>%
  ungroup() %>%
  distinct(Sample_ID, Diagnosis) %>%
  pull(Sample_ID)

# Write subjects with NaN di_gaussian round 1 data to a CSV to read into spyder
write.csv(SCZ_di_gauss_mod_NaN,
          "UCLA_Schizophrenia_di_gaussian_NaN_subjects.csv",
          row.names=F)

# Load round 1 ABIDE ASD di_gaussian data
SCZ_di_gaussian_v1 <- readRDS(paste0(SCZ_rdata_path,
                                           "UCLA_Schizophrenia_di_gaussian.Rds"))
SCZ_di_gauss_v1_NaN <- SCZ_di_gaussian_v1 %>%
  filter(SPI=="di_gaussian",
         brain_region_1 != brain_region_2) %>%
  filter(is.na(value)) %>%
  group_by(Sample_ID, Diagnosis) %>%
  summarise(num_NA = n()) %>%
  ungroup() %>%
  distinct(Sample_ID, Diagnosis) %>%
  pull(Sample_ID)

SCZ_di_gauss_v1_NaN
# This time, there were no NaNs returned for any of these subjects.

### ABIDE ASD
ASD_pyspi14_mod %>%
  filter(SPI=="di_gaussian",
         brain_region_1 != brain_region_2) %>%
  filter(is.na(value)) %>%
  group_by(Sample_ID, Diagnosis) %>%
  summarise(num_NA = n()) %>%
  ungroup() %>%
  summarise(num_full_NA = sum(num_NA == 2256),
            num_nonfull_NA = sum(num_NA < 2256))

# How many of these subjects with NaN were originally NaN?
ASD_di_gauss_og_NaN <- ASD_pyspi14 %>%
  filter(SPI=="di_gaussian",
         brain_region_1 != brain_region_2) %>%
  filter(is.na(value)) %>%
  group_by(Sample_ID, Diagnosis) %>%
  summarise(num_NA = n()) %>%
  ungroup() %>%
  distinct(Sample_ID, Diagnosis) %>%
  pull(Sample_ID)

ASD_di_gauss_mod_NaN <- ASD_pyspi14_mod %>%
  filter(SPI=="di_gaussian",
         brain_region_1 != brain_region_2) %>%
  filter(is.na(value)) %>%
  group_by(Sample_ID, Diagnosis) %>%
  summarise(num_NA = n()) %>%
  ungroup() %>%
  distinct(Sample_ID, Diagnosis) %>%
  pull(Sample_ID)

# Find number of subjects that yielded NaN originally but no longer after
# di_gaussian was run first
length(ASD_di_gauss_og_NaN[!(ASD_di_gauss_og_NaN %in% ASD_di_gauss_mod_NaN)])
# Find number of subjects that now yield NaN that didn't yield NaN previously
length(ASD_di_gauss_mod_NaN[!(ASD_di_gauss_mod_NaN %in% ASD_di_gauss_og_NaN)])

# Write subjects with NaN di_gaussian round 1 data to a CSV
write.csv(ASD_di_gauss_mod_NaN,
          "ABIDE_ASD_di_gaussian_NaN_subjects.csv",
          row.names=F)

# Load round 1 ABIDE ASD di_gaussian data
ASD_di_gaussian_v1 <- readRDS(paste0(ASD_rdata_path,
                                           "ABIDE_ASD_di_gaussian.Rds"))
ASD_di_gauss_v1_NaN <- ASD_di_gaussian_v1 %>%
  filter(SPI=="di_gaussian",
         brain_region_1 != brain_region_2) %>%
  filter(is.na(value)) %>%
  group_by(Sample_ID, Diagnosis) %>%
  summarise(num_NA = n()) %>%
  ungroup() %>%
  distinct(Sample_ID, Diagnosis) %>%
  pull(Sample_ID)

ASD_di_gauss_v1_NaN

# Write subjects with NaN di_gaussian to be run in round 2 data to a CSV
write.csv(ASD_di_gauss_v1_NaN,
          "ABIDE_ASD_di_gaussian_NaN_subjects_v2.csv",
          row.names=F)

# Load round 2 ABIDE ASD di_gaussian data
ASD_di_gaussian_v2 <- readRDS(paste0(ASD_rdata_path,
                                     "ABIDE_ASD_di_gaussian_v2.Rds"))

ASD_di_gauss_v2_NaN <- ASD_di_gaussian_v2 %>%
  filter(SPI=="di_gaussian",
         brain_region_1 != brain_region_2) %>%
  filter(is.na(value)) %>%
  group_by(Sample_ID, Diagnosis) %>%
  summarise(num_NA = n()) %>%
  ungroup() %>%
  distinct(Sample_ID, Diagnosis) %>%
  pull(Sample_ID)

ASD_di_gauss_v2_NaN


################################################################################
# Integrating corrected di_gaussian data with pyspi14_mod data
# di_gaussian for UCLA SCZ, di_gaussian_v2 for ABIDE ASD

# UCLA SCZ
if (!file.exists(paste0(SCZ_rdata_path, "UCLA_Schizophrenia_pyspi14_corrected_filtered.Rds"))) {
  # Find subjects in SCZ_di_gaussian_v1
  SCZ_corrected_subjects <- unique(SCZ_di_gaussian_v1$Sample_ID)
  # Replace original di_gaussian values with corrected (non-NaN) ones
  SCZ_pyspi14_mod_corr <- SCZ_pyspi14_mod %>%
    filter(!(Sample_ID %in% SCZ_corrected_subjects & SPI=="di_gaussian")) %>%
    plyr::rbind.fill(., SCZ_di_gaussian_v1)
  # Save to corrected Rds file
  saveRDS(SCZ_pyspi14_mod_corr, 
          paste0(SCZ_rdata_path, "UCLA_Schizophrenia_pyspi14_corrected_filtered.Rds"))
  
  # Save sample data post-filtering to an `.Rds` file:
  filtered_sample_info <- SCZ_pyspi14_mod_corr %>%
    distinct(Sample_ID)
  saveRDS(filtered_sample_info, file=paste0(SCZ_rdata_path,
                                            "UCLA_Schizophrenia_filtered_sample_info_pyspi14_corrected.Rds"))
  
  # z-score the corrected data
  SCZ_pyspi14_mod_corr <- SCZ_pyspi14_mod_corr %>%
    dplyr::rename("names"="SPI", "values"="value")
  
  SCZ_pyspi14_mod_corr_z <- normalise_feature_frame(SCZ_pyspi14_mod_corr, 
                                                    names_var = "names",
                                                    values_var = "values", 
                                                    method = "z-score")
  # Save the z-scored correct data to Rds file
  saveRDS(SCZ_pyspi14_mod_corr_z, 
          paste0(SCZ_rdata_path, "UCLA_Schizophrenia_pyspi14_corrected_filtered_zscored.Rds"))
}

# ABIDE ASD
if (!file.exists(paste0(ASD_rdata_path, "ABIDE_ASD_pyspi14_corrected_filtered.Rds"))) {
  # Find subjects in ASD_di_gaussian_v1
  ASD_corrected_subjects_v1 <- unique(ASD_di_gaussian_v1$Sample_ID)
  ASD_corrected_subjects_v2 <- unique(ASD_di_gaussian_v2$Sample_ID)
  
  # Replace original di_gaussian values with corrected (non-NaN) ones
  ASD_pyspi14_mod_corr <- ASD_pyspi14_mod %>%
    filter(!(Sample_ID %in% ASD_corrected_subjects_v1 & SPI=="di_gaussian")) %>%
    plyr::rbind.fill(., ASD_di_gaussian_v1) %>%
    # Iteratively replace the 5 subjects for whom di_gaussian still returned
    # NaNs from the first run with only di_gaussian in spyder
    filter(!(Sample_ID %in% ASD_corrected_subjects_v2 & SPI=="di_gaussian")) %>%
    plyr::rbind.fill(., ASD_di_gaussian_v2)
  # Save to corrected Rds file
  saveRDS(ASD_pyspi14_mod_corr, 
          paste0(ASD_rdata_path, "ABIDE_ASD_pyspi14_corrected_filtered.Rds"))
  
  # Save sample data post-filtering to an `.Rds` file:
  filtered_sample_info <- ASD_pyspi14_mod_corr %>%
    distinct(Sample_ID)
  saveRDS(filtered_sample_info, file=paste0(ASD_rdata_path,
                                            "ABIDE_ASD_filtered_sample_info_pyspi14_corrected.Rds"))
  
  # z-score the corrected data
  ASD_pyspi14_mod_corr <- ASD_pyspi14_mod_corr %>%
    dplyr::rename("names"="SPI", "values"="value")
  
  ASD_pyspi14_mod_corr_z <- normalise_feature_frame(ASD_pyspi14_mod_corr, 
                                                    names_var = "names",
                                                    values_var = "values", 
                                                    method = "z-score")
  # Save the z-scored correct data to Rds file
  saveRDS(ASD_pyspi14_mod_corr_z, 
          paste0(ASD_rdata_path, "ABIDE_ASD_pyspi14_corrected_filtered_zscored.Rds"))
  
}

################################################################################
# di_gaussian robustness analysis
# Subject: sub-10159
# brain region 1: ctx-lh-bankssts
# brain region 2: ctx-lh-entorhinal
# noise processing: AROMA+2P+GMR

# Subset data down to just those two brain regions for sub-10159
di_gaussian_robustness_TS_data <- SCZ_TS %>%
  filter(Sample_ID == "sub-10159",
         Noise_Proc == "AROMA+2P+GMR",
         Brain_Region %in% c("ctx-lh-bankssts", "ctx-lh-entorhinal")) %>%
  dplyr::select(Brain_Region, timepoint, values) %>%
  pivot_wider(id_cols = Brain_Region,
              names_from = "timepoint",
              values_from = "values") %>%
  dplyr::select(-Brain_Region)
write.table(di_gaussian_robustness_TS_data, 
            file="sub-10159_lh_bankssts_lh_entorhinal.csv", 
            sep=",", col.names=F, row.names=F)

# Once data has been processed with pyspi for di_gaussian 1000x,
# visualise it here
SCZ_pydata_path <- paste0(SCZ_data_path, "raw_data/pydata/")
processed_di_gauss_data_1000x <- read.csv(paste0(SCZ_pydata_path,
                                                 "sub-10159_lh_bankssts_lh_entorhinal_di_gaussian.csv"),
                                          header=F)
colnames(processed_di_gauss_data_1000x) <- "di_gaussian"
processed_di_gauss_data_1000x$Iteration <- 1:nrow(processed_di_gauss_data_1000x)

processed_di_gauss_data_1000x %>%
  ggplot(data=., mapping=aes(x=di_gaussian)) +
  geom_histogram(fill="turquoise3") +
  ggtitle("di_gaussian for sub-10159\nleft bankssts --> entorhinal") +
  ylab("# Iterations") +
  xlab("di_gaussian value") +
  theme(plot.title = element_text(hjust=0.5))
ggsave(paste0(plot_path, "di_gaussian_robustness_1000x_histogram.png"), 
       width=6, height=4, units="in", dpi=300, bg="white")