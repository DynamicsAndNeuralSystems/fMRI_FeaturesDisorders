################################################################################
# Define study/data paths
################################################################################
library(feather)
library(tidyverse)
library(glue)
library(icesTAF)
library(cowplot)
library(RNifti)
library(oro.nifti)
library(neurobase)
theme_set(theme_cowplot())

python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"
reticulate::use_python(python_to_use)
library(reticulate)

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

github_dir <- "~/Library/CloudStorage/OneDrive-TheUniversityofSydney(Students)/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/univariate_results/")
TAF::mkdir(plot_path)

# Define atlases to read in
UCLA_CNP_atlas <- "~/data/neuroimaging_atlases/mni152_space/modded/aparc+aseg.nii.gz"
ABIDE_atlas <- "~/data/neuroimaging_atlases/mni152_space/original/ho_roi_atlas.nii.gz"
aparc_aseg_lookup <- read.table("~/data/neuroimaging_atlases/FreeSurferLUT.txt", header=T)

UCLA_CNP_atlas <- pyarrow_feather$read_feather("~/data/neuroimaging_atlases/mni152_space/modded/aparc+aseg_atlas.feather") %>%
  dplyr::select(Voxel_x, Voxel_y, Voxel_z, Value)%>%
  left_join(., aparc_aseg_lookup) %>%
  dplyr::select(-Value)

UCLA_CNP_volumes <- UCLA_CNP_atlas %>%
  group_by(Brain_Region) %>%
  summarise(Num_Voxels = n(),
            Volume_mm3 = 0.737463*n())

# Load univariate classification results across all folds
univariate_balanced_accuracy_AUC_all_folds <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_balanced_accuracy_AUC_all_folds.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set, kernel==SVM_kernel)
# Compute mean + SD performance across all folds
univariate_balanced_accuracy <- univariate_balanced_accuracy_AUC_all_folds %>%
  group_by(Study, Comparison_Group, Univariate_Feature_Set, Analysis_Type, group_var, kernel) %>%
  reframe(Balanced_Accuracy_Across_Folds = mean(Balanced_Accuracy, na.rm=T),
          Balanced_Accuracy_Across_Folds_SD = sd(Balanced_Accuracy, na.rm=T),
          ROC_AUC_Across_Folds = mean(ROC_AUC, na.rm=T),
          ROC_AUC_Across_Folds_SD = sd(ROC_AUC, na.rm=T))
# Load p-values
univariate_p_values <- pyarrow_feather$read_feather(glue("{data_path}/UCLA_CNP_ABIDE_ASD_univariate_mixedsigmoid_scaler_empirical_p_values.feather")) %>%
  filter(Univariate_Feature_Set == univariate_feature_set) %>%
  dplyr::select(-Balanced_Accuracy_Across_Repeats, -Balanced_Accuracy_Across_Repeats_SD, 
                -ROC_AUC_Across_Repeats, -ROC_AUC_Across_Repeats_SD) %>%
  left_join(., univariate_balanced_accuracy)

univariate_p_values %>%
  filter(Analysis_Type == "Univariate_Brain_Region",
         Study == "UCLA_CNP") %>%
  dplyr::rename("Brain_Region" = "group_var") %>%
  left_join(., UCLA_CNP_volumes) %>%
  ggplot(data=., mapping=aes(x=Volume_mm3, y=Balanced_Accuracy_Across_Folds)) +
  geom_point() + 
  facet_grid(Comparison_Group ~ ., scales="free")