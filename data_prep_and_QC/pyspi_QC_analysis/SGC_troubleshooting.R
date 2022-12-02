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
# Frequency domain visualization
SGC_subject <- SCZ_TS %>%
  filter(Sample_ID == "sub-10527",
         Brain_Region %in% c("ctx-lh-rostralanteriorcingulate", 
                             "ctx-lh-caudalmiddlefrontal")) %>%
  mutate(Brain_Region = factor(Brain_Region, levels=c("ctx-lh-rostralanteriorcingulate", 
                                                        "ctx-lh-caudalmiddlefrontal"))) %>%
  group_by(Brain_Region) %>%
  summarise(Region_FFT = abs(fft(values)/sqrt(128))^2) %>%
  summarise(Power = (4/128)*Region_FFT[1:65],
            frequency = (0:64)/128)

SGC_subject %>%
  ggplot(data=., mapping=aes(x=frequency, y=Power, group=Brain_Region, color=Brain_Region)) +
  geom_vline(xintercept=0.25, linetype=2, alpha=0.7, color="black") +
  geom_line() +
  labs(color="Brain Region") +
  guides(color=guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.position = "bottom") 