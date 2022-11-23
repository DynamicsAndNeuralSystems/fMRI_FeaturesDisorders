################################################################################
# Load libraries
################################################################################

library(tidyverse)
library(icesTAF)
library(cowplot)
library(knitr)
library(kableExtra)
theme_set(theme_cowplot())

################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/QC/")
TAF::mkdir(plot_path)

UCLA_data_path <- "~/data/UCLA_Schizophrenia/"
UCLA_rdata_path <- paste0(UCLA_data_path, "processed_data/Rdata/")

# Load raw (non-normalized) pyspi14 data for UCLA Schizophrenia
UCLA_pyspi14 <- readRDS(paste0(UCLA_rdata_path, 
                          "UCLA_Schizophrenia_pyspi14_filtered.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR",
         Diagnosis %in% c("Control", "Schizophrenia"))

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
plot_SPI_na <- function(this_SPI, SPI_name, pyspi_data) {
  pyspi_data %>%
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
            SPI_name="sgc_full_range",
            pyspi_data = UCLA_pyspi14)
plot_SPI_na(this_SPI="di_gaussian",
            SPI_name="di_gaussian",
            pyspi_data = UCLA_pyspi14)

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


################################################################################

# I next tried re-running pyspi on the physics cluster with
# sgc_nonparametric_mean_fs-1_fmin-0-25_fmax-0-5 instead of sgc_nonparametric_mean_fs-1_fmin-0_fmax-0-5
# And I put di_gaussian at the top of the SPI computation list to see if
# Changing the order around made the difference

# I named this dataset with pyspi14_mod instead of pyspi14
UCLA_pyspi14_mod <- readRDS(paste0(UCLA_rdata_path, 
                               "UCLA_Schizophrenia_pyspi14_mod_filtered.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR")

# Check distribution of remaining NA for SGC
UCLA_pyspi14_mod %>%
  filter(SPI=="sgc_nonparametric_mean_fs-1_fmin-0-25_fmax-0-5",
         brain_region_1 != brain_region_2) %>%
  filter(is.na(value)) %>%
  arrange(brain_region_1, brain_region_2, Sample_ID) %>%
  dplyr::select(-SPI, -Noise_Proc) %>%
  kable(.) %>%
  kable_styling(full_width=F)

# Find # NaN for SGC for each brain region
UCLA_pyspi14_mod %>%
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

# di_gaussian redo
# There are some subjects with just a few regions returning NaN for di_gaussian,
# which isn't overly concerning
UCLA_pyspi14_mod %>%
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
            pyspi_data = UCLA_pyspi14_mod)

UCLA_pyspi14_mod %>%
  filter(SPI=="di_gaussian",
         brain_region_1 != brain_region_2) %>%
  filter(is.na(value)) %>%
  group_by(Sample_ID, Diagnosis) %>%
  summarise(num_NA = n()) %>%
  ungroup() %>%
  pull(Sample_ID)


# I re-ran just di_gaussian for these 13 subjects with one or more NaN directly 
# within spyder on the physics cluster
UCLA_pyspi14_mod_di_gaussian <- readRDS(paste0(UCLA_rdata_path, 
                                   "UCLA_Schizophrenia_pyspi14_mod_di_gaussian_filtered.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR")

UCLA_pyspi14_mod_di_gaussian %>%
  filter(SPI=="di_gaussian",
         brain_region_1 != brain_region_2) %>%
  filter(is.na(value)) %>%
  group_by(Sample_ID, Diagnosis) %>%
  summarise(num_NA = n()) %>%
  ungroup() %>%
  arrange(desc(num_NA)) %>%
  kable(.) %>%
  kable_styling(full_width=F)

# This time, there were no NaNs returned for any of these subjects.