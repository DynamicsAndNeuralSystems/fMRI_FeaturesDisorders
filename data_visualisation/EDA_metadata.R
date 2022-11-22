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
source(paste0(github_dir, "helper_functions/classification/Linear_SVM.R"))
plot_path <- paste0(github_dir, "plots/EDA/")
icesTAF::mkdir(plot_path)

SCZ_data_path <- "~/data/UCLA_Schizophrenia/"
SCZ_rdata_path <- paste0(SCZ_data_path, "processed_data/Rdata/")
ASD_data_path <- "~/data/ABIDE_ASD/"
ASD_rdata_path <- paste0(ASD_data_path, "processed_data/Rdata/")

# Load subject metadata
SCZ_subject_metadata <- readRDS(paste0(SCZ_data_path, "UCLA_Schizophrenia_sample_metadata.Rds")) %>%
  dplyr::filter(Diagnosis %in% c("Control", "Schizophrenia"))
ASD_subject_metadata <- readRDS(paste0(ASD_data_path, "ABIDE_ASD_sample_metadata.Rds"))

# Summary of SCZ dataset 
SCZ_subject_metadata %>%
  group_by(Diagnosis) %>%
  summarise(N = n(),
            Percent_Female = 100*sum(gender=="F", na.rm=T)/n(),
            Num_Female = sum(gender=="F", na.rm=T),
            Age_Mean = mean(age, na.rm=T),
            Age_SD = sd(age, na.rm=T))

# Summary of ASD dataset 
ASD_subject_metadata %>%
  mutate(Diagnosis = factor(Diagnosis, levels = c("Control", "ASD"))) %>%
  group_by(Diagnosis) %>%
  summarise(N = n(),
            Percent_Female = 100*sum(sex=="F", na.rm=T)/n(),
            Num_Female = sum(sex=="F", na.rm=T),
            Age_Mean = mean(age, na.rm=T),
            Age_SD = sd(age, na.rm=T))

# Visualise age and sex distribution by diagnosis
SCZ_subject_metadata %>%
  mutate(Diagnosis = factor(Diagnosis, levels = c("Control", "Schizophrenia"))) %>%
  ggplot(data = ., mapping=aes(x = gender, y = age, fill = Diagnosis)) +
  geom_violin(position = position_dodge(0.9)) +
  geom_boxplot(color = "black", width=0.1, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#00B06D", "#737373")) +
  ggtitle("Schizophrenia Age\nvs. Sex and Diagnosis") +
  ylab("Age") +
  xlab("Sex") +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust=0.5))
ggsave(paste0(plot_path, "SCZ_Age_vs_Sex_and_Diagnosis.png"),
       width=5, height=4, units="in", dpi=300,
       bg = "white")

ASD_subject_metadata %>%
  mutate(Diagnosis = factor(Diagnosis, levels = c("Control", "ASD"))) %>%
  ggplot(data = ., mapping=aes(x = sex, y = age, fill = Diagnosis)) +
  geom_violin(position = position_dodge(0.9)) +
  geom_boxplot(color = "black", width=0.1, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#00B06D", "#737373")) +
  ggtitle("ASD Age\nvs. Sex and Diagnosis") +
  ylab("Age") +
  xlab("Sex") +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust=0.5))
ggsave(paste0(plot_path, "ASD_Age_vs_Sex_and_Diagnosis.png"),
       width=5, height=4, units="in", dpi=300,
       bg = "white")


# Look into sites for ASD data
ASD_subject_metadata %>%
  mutate(site = ifelse(site > 20, site - 20, site)) %>%
  group_by(Diagnosis, site) %>%
  count() %>%
  pivot_wider(id_cols = site,
              names_from = Diagnosis,
              values_from = n) %>%
  mutate(ASD = ifelse(is.na(ASD), 0, ASD),
         Control = ifelse(is.na(Control), 0, Control))


SVM_res <- readRDS(paste0(ASD_rdata_path, 
                          "ROI_wise_CV_linear_SVM_catch22_inv_prob.Rds")) 

SVM_res %>%
  distinct(Sample_ID, Actual_Diagnosis) %>%
  group_by(Actual_Diagnosis) %>%
  count()