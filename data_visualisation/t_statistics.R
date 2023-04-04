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
library(knitr)
library(kableExtra)
library(patchwork)
library(broom)
library(colorspace)
library(see)
library(ggridges)
library(scales)
library(ggseg)
library(ggsegHO)
library(igraph)
library(ggraph)
library(LaCroixColoR)
theme_set(theme_cowplot())

# Source visualisation script
source(glue("{github_dir}/data_visualisation/Manuscript_Draft_Visualisations_Helper.R"))

# Load in univariate time-series feature info
catch24_info <- read.csv(glue("{github_dir}/data_visualisation/catch24_info.csv"))
pyspi14_info <- read.csv(glue("{github_dir}/data_visualisation/SPI_info.csv"))

# Load study metadata
UCLA_CNP_metadata <- pyarrow_feather$read_feather("~/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather") 
ABIDE_ASD_metadata <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_sample_metadata.feather") 

# Load brain region info
UCLA_CNP_brain_region_info <- read.csv("~/data/UCLA_CNP/study_metadata/UCLA_CNP_Brain_Region_info.csv")
ABIDE_ASD_brain_region_info <- read.csv("~/data/ABIDE_ASD/study_metadata/ABIDE_ASD_Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv")
region_node_to_from <- read.csv("~/data/TS_feature_manuscript/node_to_from_structure.csv")

# Load SPI info

################################################################################
# Hyper vs hypo connectivity maps



scz_example_data <- t_stats_pyspi14_whole_brain %>%
  ungroup() %>%
  group_by(Region_Pair, Comparison_Group) %>%
  summarise(t_stat_sum = mean(statistic)) %>%
  arrange(t_stat_sum) %>%
  separate(Region_Pair, into=c("from", "to"), sep="_") %>%
  filter(Comparison_Group=="SCZ") %>%
  sample_n(100)

# Edges are defined as cortical lobe --> specific ROI connection
edges <- region_node_to_from %>%
  filter(Study == "UCLA_CNP") %>% 
  distinct() %>%
  dplyr::select(-Study)

# ROIs don't include the origin --> cortical lobe connection
rois <- edges %>% filter(!(to %in% c("Cingulate", "Frontal", "Insula",
                                     "Occipital", "Parietal", "Temporal", "Subcortex")))

# Create a dataframe of vertices, one line per object in the ROI cortical lobe hierarchy
vertices = data.frame(name = unique(c(as.character(edges$from), as.character(edges$to))))
vertices$group <- edges$from[match(vertices$name, edges$to)]

# Create an igraph object
mygraph <- graph_from_data_frame(d=edges, vertices=vertices)

# connect = dataframe of pairwise correlations between cortical ROIs
connect <- scz_example_data %>%
  rename("value" = "t_stat_sum") %>%
  arrange(from, to)

# mygraph = igraph object linking each cortical ROI
# convert to a circular dendrogram-shaped ggraph object
p <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  theme_void()

from <- match(connect$from, vertices$name)
to <- match(connect$to, vertices$name)

p <- p +  geom_conn_bundle(data = get_con(from = from, to = to, 
                                          corr=connect$value), 
                           tension=0.5, width=3,
                           aes(alpha=corr), color="lightsteelblue3")  +
  labs(edge_alpha="Pearson\nCorrelation")

# Add leaf nodes showing the cortical ROI, color by cortical lobe
p <- p + geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group),   
                         size=3) +
  scale_color_manual(values=c(lacroix_palette("PassionFruit"), "red")) +
  labs(color="Cortex") +
  # geom_node_text(aes(x = x*1.2, y=y*1.2, filter = leaf, label=name,
  #                    color=group))+
  theme_void() + 
  theme(plot.title=element_text(size=14, face="bold", hjust=0.5),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(10,10,10,10))

p