################################################################################
# Load libraries
################################################################################

library(tidyverse)
library(icesTAF)
library(cowplot)
theme_set(theme_cowplot())

################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/QC/")
icesTAF::mkdir(plot_path)

UCLA_data_path <- "~/data/UCLA_Schizophrenia/"

# Load raw (non-normalized) pyspi14 data for UCLA Schizophrenia
UCLA_pyspi14 <- readRDS(paste0(UCLA_rdata_path, 
                          "UCLA_Schizophrenia_pyspi14_filtered.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR")

# Find SPIs that yielded NaN for at least one region pair and subject
UCLA_pyspi14 %>%
  filter(brain_region_1 != brain_region_2) %>%
  filter(is.na(value)) %>%
  distinct(SPI) %>%
  pull(SPI)

# Two SPIs yielded NaN for at least one region pair and subject:
# sgc_nonparametric_mean_fs-1_fmin-0_fmax-0-5
# di_gaussian

# Heatmap for each of these showing the subjects (rows)  and brain region pairs (columns)
plot_SPI_na <- function(this_SPI, SPI_name) {
  UCLA_pyspi14 %>%
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
  ggsave(paste0(plot_path, SPI_name, "_NaN_heatmap_UCLA_Schizophrenia.png"),
         width=12, height=3, units="in", dpi=1200, bg="white")
}
plot_SPI_na(this_SPI="sgc_nonparametric_mean_fs-1_fmin-0_fmax-0-5",
            SPI_name="sgc_full_range")
plot_SPI_na(this_SPI="di_gaussian",
            SPI_name="di_gaussian")

# Load raw time-series data for UCLA Schizophrenia
UCLA_TS <- readRDS(paste0(UCLA_data_path, 
                          "raw_data/UCLA_Schizophrenia_fMRI_TS.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR")

# First pick a region pair that yielded NaN for sgc full range from pyspi
# Focusing on sub-10527 with AROMA+2P+GMR noise processing
UCLA_TS %>%
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
UCLA_TS %>%
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