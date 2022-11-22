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

# Load raw time-series data for UCLA Schizophrenia
UCLA_TS <- readRDS(paste0(UCLA_data_path, 
                          "raw_data/UCLA_Schizophrenia_fMRI_TS.Rds")) %>%
  filter(Noise_Proc == "AROMA+2P+GMR")

# First pick a region pair that yielded NaN from pyspi
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

UCLA_pyspi14 %>%
  filter(Sample_ID == "sub-10527",
         names == "sgc_nonparametric_mean_fs-1_fmin-0_fmax-0-5",
         brain_region_1 == "ctx-lh-caudalmiddlefrontal",
         brain_region_2 == "ctx-lh-rostralanteriorcingulate") %>%
  pull(values)


# read in temp.csv
temp_res = read.csv("temp.csv") %>%
  dplyr::select(-X)

temp_res %>%
  mutate(brain_region_1 = 1 + as.numeric(str_replace_all(brain_region_1, "proc-", "")),
         brain_region_2 = 1 + as.numeric(str_replace_all(brain_region_2, "proc-", ""))) %>%
  filter(brain_region_1 != brain_region_2) %>%
  filter(is.na(value)) %>%
  group_by(SPI) %>%
  count()