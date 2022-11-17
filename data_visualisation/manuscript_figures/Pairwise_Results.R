

################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/")
icesTAF::mkdir(plot_path)

SCZ_data_path <- "~/data/UCLA_Schizophrenia/"
SCZ_rdata_path <- paste0(SCZ_data_path, "processed_data/Rdata/")
ASD_data_path <- "~/data/ABIDE_ASD/"
ASD_rdata_path <- paste0(ASD_data_path, "processed_data/Rdata/")

raw_ASD_pyspi_res <- readRDS(paste0(ASD_rdata_path, "ABIDE_ASD_pyspi14_filtered.Rds"))

# Get number of brain region pairs
raw_ASD_pyspi_res %>%
  filter(brain_region_1 != brain_region_2) %>%
  distinct(brain_region_1, brain_region_2) %>%
  nrow()

sgc_res <- raw_ASD_pyspi_res %>%
  filter(brain_region_1 != brain_region_2) %>%
  filter(SPI=="sgc_nonparametric_mean_fs-1_fmin-0_fmax-0-5") 

sgc_res_NA <-  sgc_res %>%
  filter(is.na(value))

sgc_res_NA %>%
  distinct(brain_region_1, brain_region_2)

sgc_res_NA %>%
  mutate(region_pair = paste0(brain_region_1, "_", brain_region_2)) %>%
  ggplot(data=., mapping=aes(x=region_pair, y=Sample_ID, fill=value)) +
  geom_tile() +
  theme(axis.text = element_blank(),
        legend.position = "none")