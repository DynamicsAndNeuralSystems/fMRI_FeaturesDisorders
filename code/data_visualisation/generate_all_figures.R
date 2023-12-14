
################################################################################
# Define study/data paths
################################################################################

# Set desired conda environment version of python to use pyarrow.feather
python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3" # TOCHANGE
reticulate::use_python(python_to_use)
library(reticulate)
pyarrow_feather <- import("pyarrow.feather")

# Define where data is stored
data_path <- "~/data/TS_feature_manuscript" # TOCHANGE

pairwise_feature_set <- "pyspi14"
univariate_feature_set <- "catch25"
SVM_kernel <- "Linear" # Change if you used a different kernel
study_group_df <- data.frame(Study = c(rep("UCLA_CNP", 3), "ABIDE_ASD"),
                             Noise_Proc = c(rep("AROMA+2P+GMR",3), "FC1000"),
                             Comparison_Group = c("Schizophrenia", "Bipolar", "ADHD", "ASD"),
                             Group_Nickname = c("SCZ", "BP", "ADHD", "ASD"))

################################################################################
# Load libraries
################################################################################

library(broom)
library(circlize)
library(colorspace)
library(ComplexHeatmap)
library(correctR)
library(cowplot)
library(dendextend)
library(factoextra)
library(FactoMineR)
library(feather)
library(ggnewscale)
library(ggpp)
library(ggpubr)
library(ggraph)
library(ggridges)
library(ggseg)
library(ggsegDefaultExtra)
library(ggsegHO)
library(ggsignif)
library(glue)
library(icesTAF)
library(igraph)
library(LaCroixColoR)
library(patchwork)
library(R.matlab)
library(RColorBrewer)
library(reticulate)
library(scales)
library(see)
library(splitstackshape)
library(theft)
library(tidyverse)

theme_set(theme_cowplot())

################################################################################
# Figure 1: Methods
################################################################################