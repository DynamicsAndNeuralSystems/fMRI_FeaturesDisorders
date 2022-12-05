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

SCZ_data_path <- "~/data/UCLA_CNP/"
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
                         "raw_data/UCLA_CNP_fMRI_TS.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR")
# Load raw time-series data for ABIDE ASD
ASD_TS <- readRDS(paste0(ASD_data_path, 
                         "raw_data/ABIDE_ASD_fMRI_TS.Rds")) %>%
  filter(Noise_Proc == "FC1000")

# Load raw (non-normalized) SGC data for UCLA Schizophrenia
SCZ_SGC <- readRDS(paste0(SCZ_rdata_path, 
                              "UCLA_Schizophrenia_pyspi14_SGC_filtered.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR",
         brain_region_from != brain_region_to,
         Diagnosis %in% c("Control", "Schizophrenia"))
ASD_pyspi14 <- readRDS(paste0(ASD_rdata_path, 
                              "ABIDE_ASD_pyspi14_filtered.Rds")) %>%
  filter(Noise_Proc == "FC1000",
         brain_region_from != brain_region_to,
         Diagnosis %in% c("Control", "ASD"))

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



################################################################################
# Find NaNs in data
################################################################################

# Find an example that is NaN for all 3 SGC ranges (full, lower, upper)
SGC_NaN_for_all_freq_ranges <- SCZ_SGC %>%
  filter(is.na(value)) %>%
  group_by(Sample_ID, brain_region_from, brain_region_to) %>%
  count() %>%
  filter(n==3) %>%
  dplyr::select(-n) 

SGC_NaN_for_all_freq_ranges %>%
  kable() %>%
  kable_styling(full_width=F)

# Focus on sub-10274 ctx-lh-lingual --> ctx-rh-lingual
SCZ_TS %>%
  filter(Sample_ID == "sub-10274",
         Brain_Region %in% c("ctx-lh-lingual",
                             "ctx-rh-lingual")) %>%
  mutate(Brain_Region = factor(Brain_Region,
                               levels = c("ctx-lh-lingual",
                                          "ctx-rh-lingual"))) %>%
  ggplot(data=., mapping=aes(x=timepoint, y=values, color=Brain_Region)) +
  geom_line(size=1, alpha=0.7) +
  scale_color_manual(values=c("darkgreen", "darkorchid")) +
  ggtitle("sub-10274 with NaN for\nAll SGC frequencies") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust=0.5)) +
  labs(color="Brain Region") +
  xlab("Time point") +
  ylab("AROMA+2P+GMR Signal Value") + 
  guides(color=guide_legend(nrow=2,byrow=TRUE))
ggsave(paste0(plot_path, "sgc_with_NaN_all_frequencies_sub-10274_lingual_TS.png"),
       width=6, height=4, units="in", dpi=300, bg="white")

# Plot this data in the frequency domain
SGC_full_NaN <- SCZ_TS %>%
  filter(Sample_ID == "sub-10274",
         Brain_Region %in% c("ctx-lh-lingual",
                             "ctx-rh-lingual")) %>%
  mutate(Brain_Region = factor(Brain_Region,
                               levels = c("ctx-lh-lingual",
                                          "ctx-rh-lingual"))) %>%
  group_by(Brain_Region) %>%
  summarise(Region_FFT = abs(fft(values)/sqrt(128))^2) %>%
  summarise(Power = (4/128)*Region_FFT[1:65],
            frequency = (0:64)/128)

SGC_full_NaN %>%
  ggplot(data=., mapping=aes(x=frequency, y=Power, group=Brain_Region, color=Brain_Region)) +
  geom_vline(xintercept=0.25, linetype=2, alpha=0.7, color="black") +
  geom_line(size=1, alpha=0.7) +
  scale_color_manual(values=c("darkgreen", "darkorchid")) +
  labs(color="Brain Region") +
  guides(color=guide_legend(nrow=1,byrow=TRUE)) +
  theme(legend.position = "bottom") 
ggsave(paste0(plot_path, "sgc_with_NaN_all_frequencies_sub-10274_lingual_FFT.png"),
       width=6, height=4, units="in", dpi=300, bg="white")

# I'll write this data to a CSV to play with in spyder directly

# Subject: sub-10274
# brain region 1: ctx-lh-lingual
# brain region 2: ctx-rh-lingual
# noise processing: AROMA+2P+GMR
subset_to_subject_region(input_TS_data = SCZ_TS, 
                         pydata_path = SCZ_pydata_path,
                         sample_ID = "sub-10274",
                         brain_regions = c("ctx-lh-lingual", 
                                           "ctx-rh-lingual"),
                         noise_proc = SCZ_noise_proc,
                         file_name = "sub-10274_lh_lingual_rh_lingual.csv")


# Plot periodogram with red lines to indicate where SGC was NaN from spectral-connectivity
SGC_full_NaN %>%
  ggplot(data=., mapping=aes(x=frequency, y=Power, group=Brain_Region, color=Brain_Region)) +
  geom_vline(xintercept=0.25, linetype=2, alpha=0.7, color="black") +
  geom_line(size=1, alpha=0.7) +
  scale_color_manual(values=c("darkgreen", "darkorchid")) +
  labs(color="Brain Region") +
  guides(color=guide_legend(nrow=1,byrow=TRUE)) +
  geom_vline(xintercept = c(0.01875, 0.05625, 0.06875, 0.08125, 0.25625, 0.26875),
             color = "red", alpha=0.8, size=1) +
  theme(legend.position = "bottom") 
ggsave(paste0(plot_path, "sgc_with_NaN_all_frequencies_sub-10274_lingual_FFT_w_NaN.png"),
       width=6, height=4, units="in", dpi=300, bg="white")


################################################################################
# Frequency domain visualization

# Subject with NaN at low end of range: sub-10527, left RAC --> left CMF
SGC_low_NaN <- SCZ_TS %>%
  filter(Sample_ID == "sub-10527",
         Brain_Region %in% c("ctx-lh-rostralanteriorcingulate", 
                             "ctx-lh-caudalmiddlefrontal")) %>%
  mutate(Brain_Region = factor(Brain_Region, levels=c("ctx-lh-rostralanteriorcingulate", 
                                                      "ctx-lh-caudalmiddlefrontal"))) %>%
  group_by(Brain_Region) %>%
  summarise(Region_FFT = abs(fft(values)/sqrt(128))^2) %>%
  summarise(Power = (4/128)*Region_FFT[1:65],
            frequency = (0:64)/128)

SGC_low_NaN %>%
  ggplot(data=., mapping=aes(x=frequency, y=Power, group=Brain_Region, color=Brain_Region)) +
  geom_vline(xintercept=0.25, linetype=2, alpha=0.7, color="black") +
  geom_line() +
  labs(color="Brain Region") +
  guides(color=guide_legend(nrow=1,byrow=TRUE)) +
  theme(legend.position = "bottom") 


# Subject with NaN at high end of range: sub-10206, right precentral <--> right postcentral
SGC_high_NaN <- SCZ_TS %>%
  filter(Sample_ID == "sub-10206",
         Brain_Region %in% c("ctx-rh-precentral", 
                             "ctx-rh-postcentral")) %>%
  mutate(Brain_Region = factor(Brain_Region, levels=c("ctx-rh-precentral", 
                                                      "ctx-rh-postcentral"))) %>%
  group_by(Brain_Region) %>%
  summarise(Region_FFT = abs(fft(values)/sqrt(128))^2) %>%
  summarise(Power = (4/128)*Region_FFT[1:65],
            frequency = (0:64)/128)

SGC_high_NaN %>%
  ggplot(data=., mapping=aes(x=frequency, y=Power, group=Brain_Region, color=Brain_Region)) +
  geom_vline(xintercept=0.25, linetype=2, alpha=0.7, color="black") +
  geom_line() +
  labs(color="Brain Region") +
  guides(color=guide_legend(nrow=1,byrow=TRUE)) +
  theme(legend.position = "bottom") 