################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/t_statistics/")
TAF::mkdir(plot_path)

# python_to_use <- "~/.conda/envs/pyspi/bin/python3"
python_to_use <- "/Users/abry4213/opt/anaconda3/envs/pyspi/bin/python3"
pairwise_feature_set <- "pyspi14"
univariate_feature_set <- "catch24"
data_path <- "~/data/TS_feature_manuscript"
study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Noise_Proc = c(rep("AROMA+2P+GMR",3), "FC1000"),
                             Comparison_Group = c("Schizophrenia", "Bipolar", "ADHD", "ASD"),
                             Group_Nickname = c("SCZ", "BPD", "ADHD", "ASD"))


ABIDE_ASD_brain_region_info <- read.csv("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv")

reticulate::use_python(python_to_use)

library(reticulate)

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

################################################################################
# Load libraries
################################################################################
library(feather)
library(tidyverse)
library(glue)
library(icesTAF)
library(cowplot)
library(ggseg)
library(ggsegHO)
library(knitr)
library(kableExtra)
library(patchwork)
library(ggseg)
library(broom)
library(colorspace)
library(see)
library(ggridges)
theme_set(theme_cowplot())

# Source visualisation script
source(glue("{github_dir}/data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))

# Load in univariate time-series feature info
TS_feature_info <- read.csv(glue("{github_dir}/data_visualisation/catch24_info.csv"))

# Load study metadata
UCLA_CNP_metadata <- pyarrow_feather$read_feather("~/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather") 
ABIDE_ASD_metadata <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_sample_metadata.feather") 

# Load raw feature data
UCLA_CNP_catch24 <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_catch24_filtered.feather")  %>%
  left_join(., UCLA_CNP_metadata) %>%
  mutate(label = ifelse(str_detect(Brain_Region, "ctx-"),
                        gsub("-", "_", Brain_Region),
                        as.character(Brain_Region))) %>%
  mutate(label = gsub("ctx_", "", label))
ABIDE_ASD_catch24 <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/processed_data/ABIDE_ASD_FC1000_catch24_filtered.feather")  %>%
  left_join(., ABIDE_ASD_metadata) %>%
  left_join(., ABIDE_ASD_brain_region_info) %>%
  dplyr::rename("region" = "ggseg")


################################################################################
# Ridge plot for catch24 features' T-statistics across entire brain
T_stats_for_group <- function(comparison_group, study, group_nickname){
  res <- plyr::rbind.fill(UCLA_CNP_catch24, ABIDE_ASD_catch24) %>%
    filter(Diagnosis %in% c(comparison_group, "Control"),
           Study == study) %>%
    mutate(Diagnosis = case_when(Diagnosis == "Schizophrenia" ~ "SCZ",
                                 Diagnosis == "Bipolar" ~ "BPD",
                                 T ~ Diagnosis)) %>%
    dplyr::select(Brain_Region, names, Diagnosis, values) %>%
    mutate(Diagnosis = factor(Diagnosis, levels = c(group_nickname, "Control"))) %>%
    group_by(Brain_Region, names) %>%
    nest() %>%
    mutate(
      fit = map(data, ~ t.test(values ~ Diagnosis, data = .x)),
      tidied = map(fit, tidy)
    ) %>% 
    unnest(tidied) %>%
    dplyr::select(-data, -fit) %>%
    arrange(p.value) %>%
    ungroup() %>%
    dplyr::select(Brain_Region, names, statistic) %>%
    dplyr::rename("TS_Feature" = "names") %>% 
    mutate(Comparison_Group = group_nickname,
           Study = study)
  
  return(res)
}

t_stats_whole_brain <- 1:4 %>%
  purrr::map_df(~ T_stats_for_group(comparison_group = study_group_df$Comparison_Group[.x],
                                    study = study_group_df$Study[.x],
                                    group_nickname = study_group_df$Group_Nickname[.x]))

t_stats_whole_brain %>%
  ungroup() %>%
  left_join(., TS_feature_info) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD")))%>%
  mutate(Figure_Name = fct_reorder(Figure_Name, statistic, .fun=sd)) %>%
  ggplot(data=., mapping=aes(x=statistic, y=Figure_Name, fill=Comparison_Group, color=Comparison_Group)) +
  geom_density_ridges(alpha=0.6, scale=1.1) +
  xlab("T-statistic across\nall brain regions") +
  ylab("catch24 time-series feature") +
  scale_fill_manual(values=c("Control" = "#5BB67B", 
                             "SCZ" = "#573DC7", 
                             "BPD" = "#D5492A", 
                             "ADHD" = "#0F9EA9", 
                             "ASD" = "#C47B2F")) +
  scale_color_manual(values=c("Control" = "#5BB67B", 
                             "SCZ" = "#573DC7", 
                             "BPD" = "#D5492A", 
                             "ADHD" = "#0F9EA9", 
                             "ASD" = "#C47B2F")) +
  guides(fill = guide_legend(nrow=2),
         color = guide_legend(nrow=2)) +
  theme(legend.position = "bottom",
        axis.title = element_text(size=16),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=15),
        legend.text = element_text(size=16),
        legend.title = element_blank())
ggsave(glue("{plot_path}/catch24_feature_t_statistics_across_brain.png"),
       width=6, height=9, units="in", dpi=300)