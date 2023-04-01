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
pyspi14_info <- read.csv(glue("{github_dir}/data_visualisation/SPI_info.csv"))

################################################################################
# Ridge plot for catch24 features' T-statistics across entire brain
T_stats_for_group <- function(comparison_group, input_data, study, group_nickname){
  res <- input_data %>%
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

if (!file.exists(glue("{data_path}/univariate_catch24_t_statistics_by_brain_region.feather"))) {
  # Load raw feature data
  UCLA_CNP_catch24 <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_catch24_filtered.feather")  %>%
    left_join(., UCLA_CNP_metadata) %>%
    mutate(Study = "UCLA_CNP")
  ABIDE_ASD_catch24 <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/processed_data/ABIDE_ASD_FC1000_catch24_filtered.feather")  %>%
    left_join(., ABIDE_ASD_metadata) %>%
    left_join(., ABIDE_ASD_brain_region_info)
  
  combined_univariate_data <- plyr::rbind.fill(UCLA_CNP_catch24, ABIDE_ASD_catch24)
  
  t_stats_catch24_whole_brain <- 1:4 %>%
    purrr::map_df(~ T_stats_for_group(input_data = combined_univariate_data,
                                      comparison_group = study_group_df$Comparison_Group[.x],
                                      study = study_group_df$Study[.x],
                                      group_nickname = study_group_df$Group_Nickname[.x]))
  feather::write_feather(t_stats_catch24_whole_brain, glue("{data_path}/univariate_catch24_t_statistics_by_brain_region.feather"))
} else {
  t_stats_catch24_whole_brain <- feather::read_feather(glue("{data_path}/univariate_catch24_t_statistics_by_brain_region.feather"))
}

t_stats_catch24_whole_brain %>%
  ungroup() %>%
  left_join(., catch24_info) %>%
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
  guides(fill = guide_legend(nrow=2, byrow=T),
         color = guide_legend(nrow=2, byrow=T)) +
  scale_y_discrete(labels = wrap_format(28)) +
  theme(legend.position = "bottom",
        axis.title = element_text(size=16),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=15),
        legend.text = element_text(size=16),
        legend.title = element_blank())
ggsave(glue("{plot_path}/catch24_feature_t_statistics_across_brain.png"),
       width=5.5, height=10, units="in", dpi=300)

# Pairwise pyspi14 T-statistics
T_stats_for_group_pairwise <- function(comparison_group, input_data, study, group_nickname){
  res <- input_data %>%
    filter(Diagnosis %in% c(comparison_group, "Control"),
           Study == study) %>%
    mutate(Diagnosis = case_when(Diagnosis == "Schizophrenia" ~ "SCZ",
                                 Diagnosis == "Bipolar" ~ "BPD",
                                 T ~ Diagnosis)) %>%
    rowwise() %>%
    mutate(Region_Pair = paste0(brain_region_from, "_", brain_region_to)) %>%
    dplyr::select(Region_Pair, SPI, Diagnosis, value) %>%
    mutate(Diagnosis = factor(Diagnosis, levels = c(group_nickname, "Control"))) %>%
    group_by(Region_Pair, SPI) %>%
    nest() %>%
    mutate(
      fit = map(data, ~ t.test(value ~ Diagnosis, data = .x)),
      tidied = map(fit, tidy)
    ) %>% 
    unnest(tidied) %>%
    dplyr::select(-data, -fit) %>%
    ungroup() %>%
    dplyr::select(Region_Pair, SPI, statistic) %>%
    mutate(Comparison_Group = group_nickname,
           Study = study)
  
  return(res)
}

if (!file.exists(glue("{data_path}/pairwise_pyspi14_t_statistics_by_region_pair.feather"))) {
  UCLA_CNP_pyspi14 <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_pyspi14_filtered.feather")  %>%
    left_join(., UCLA_CNP_metadata) %>%
    filter(!is.na(Diagnosis)) %>%
    mutate(Study = "UCLA_CNP")
  ABIDE_ASD_pyspi14 <- pyarrow_feather$read_feather("~/data/ABIDE_ASD/processed_data/ABIDE_ASD_FC1000_pyspi14_filtered.feather")  %>%
    left_join(., ABIDE_ASD_metadata) %>%
    filter(!is.na(Diagnosis)) %>%
    mutate(Study = "ABIDE_ASD")
  
  t_stats_pyspi14_whole_brain <- 1:4 %>%
  purrr::map_df(~ T_stats_for_group_pairwise(input_data = plyr::rbind.fill(UCLA_CNP_pyspi14,
                                                                  ABIDE_ASD_pyspi14),
                                    comparison_group = study_group_df$Comparison_Group[.x],
                                    study = study_group_df$Study[.x],
                                    group_nickname = study_group_df$Group_Nickname[.x]))

    feather::write_feather(t_stats_pyspi14_whole_brain, glue("{data_path}/pairwise_pyspi14_t_statistics_by_region_pair.feather"))
} else {
  t_stats_pyspi14_whole_brain <- feather::read_feather(glue("{data_path}/pairwise_pyspi14_t_statistics_by_region_pair.feather"))
}

t_stats_pyspi14_whole_brain %>%
  ungroup() %>%
  left_join(., pyspi14_info) %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD")))%>%
  mutate(Nickname = fct_reorder(Nickname, statistic, .fun=sd)) %>%
  ggplot(data=., mapping=aes(x=statistic, y=Nickname, fill=Comparison_Group, color=Comparison_Group)) +
  geom_density_ridges(alpha=0.6, scale=1.1) +
  xlab("T-statistic across\nall brain regions") +
  ylab("pyspi14 time-series feature") +
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
  scale_y_discrete(labels = wrap_format(28)) +
  theme(legend.position = "bottom",
        axis.title = element_text(size=16),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=15),
        legend.text = element_text(size=16),
        legend.title = element_blank())
ggsave(glue("{plot_path}/pyspi14_feature_t_statistics_across_brain.png"),
       width=5.5, height=10, units="in", dpi=300)

# Find regions most disrupted across all pairwise connections
pairwise_t_stats_by_region_from <- t_stats_pyspi14_whole_brain %>%
  separate(Region_Pair, c("region_from", "region_to"),
           sep="_") %>%
  group_by(region_from, SPI, Study, Comparison_Group) %>%
  summarise(mean_T_magnitude = mean(abs(statistic))) %>%
  dplyr::rename("Brain_Region" = "region_from") %>%
  mutate(Direction = "from")

pairwise_t_stats_by_region_to <- t_stats_pyspi14_whole_brain %>%
  separate(Region_Pair, c("region_from", "region_to"),
           sep="_") %>%
  group_by(region_to, SPI, Study, Comparison_Group) %>%
  summarise(mean_T_magnitude = mean(abs(statistic)))%>%
  dplyr::rename("Brain_Region" = "region_to") %>%
  mutate(Direction = "to")


################################################################################
# Function to plot data in the brain
group_colors <- c("#573DC7", "#D5492A", "#0F9EA9","#C47B2F")
plot_data_per_group <- function(SPI_nickname, input_data, group_colors, min_max_df=NULL) { 
  ggseg_plot_list <- list()
  legend_list <- list()
  group_colors <- c("#573DC7", "#D5492A", "#0F9EA9","#C47B2F")
  
  for (i in 1:nrow(study_group_df)) {
    dataset_ID <- study_group_df$Study[i]
    comparison_group <- study_group_df$Group_Nickname[i]
    group_color <- group_colors[i]
    
    # Define atlas by study
    atlas <- ifelse(dataset_ID == "UCLA_CNP", "dk", "hoCort")
    
    T_data_to_plot <- input_data %>%
      filter(Study == dataset_ID,
             Comparison_Group == comparison_group) %>%
      select(where(function(x) any(!is.na(x)))) %>%
      arrange(desc(mean_T_magnitude)) %>%
      ungroup() 
    
    if (!is.null(min_max_df)) {
      min_fill = min_max_df %>% 
        filter(Study == dataset_ID,
               Comparison_Group == comparison_group) %>%
        pull(min_fill)
      max_fill = min_max_df %>% 
        filter(Study == dataset_ID,
               Comparison_Group == comparison_group) %>%
        pull(max_fill)
    } else {
      min_fill = NULL
      max_fill = NULL
    }
    
    # Plot T stat data in cortex
    dataset_ggseg <- plot_data_with_ggseg(dataset_ID=dataset_ID,
                                          atlas_name=atlas,
                                          atlas_data=get(atlas),
                                          data_to_plot=T_data_to_plot,
                                          min_fill = min_fill,
                                          max_fill = max_fill,
                                          line_color = "gray30",
                                          fill_variable="mean_T_magnitude",
                                          fill_colors = c("white", group_color)) 
    
    # Add subcortical data for UCLA CNP
    if (dataset_ID == "UCLA_CNP") {
      dataset_ggseg_subctx <- plot_data_with_ggseg(dataset_ID = dataset_ID,
                                                   atlas_name = "aseg",
                                                   atlas_data = aseg,
                                                   data_to_plot=T_data_to_plot,
                                                   min_fill = min_fill,
                                                   max_fill = max_fill,
                                                   line_color = "gray30",
                                                   fill_variable="mean_T_magnitude",
                                                   fill_colors = c("white", group_color)) +
        labs(fill = glue("Mean Absolute T-Statistic\nfor {SPI_nickname}"))
      
      # Extract just legend
      subctx_legend <- ggpubr::as_ggplot(ggpubr::get_legend(dataset_ggseg_subctx + 
                                                              theme(legend.position = "bottom",
                                                                    legend.text = element_text(size=12)) +
                                                              labs(fill = glue("Mean Absolute T-Statistic\nfor {SPI_nickname}")) +
                                                              guides(fill = guide_colorbar(title.position = "top", 
                                                                                           nrow = 1,
                                                                                           barwidth = 12, 
                                                                                           barheight = 0.75,
                                                                                           title.hjust = 0.5,
                                                                                           label.position = "bottom"))))
      
      # Extract just brain
      dataset_ggseg_subctx <- dataset_ggseg_subctx + 
        theme(legend.position = "none")
      
      # Append to list
      dataset_ggseg <- dataset_ggseg + theme(legend.position = "none")
      ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg)
      ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg_subctx)
      legend_list <- list.append(legend_list, subctx_legend)
    } else {
      # For ABIDE
      # Extract just legend
      ctx_legend <- ggpubr::as_ggplot(ggpubr::get_legend(dataset_ggseg + 
                                                           theme(legend.position = "bottom",
                                                                 legend.text = element_text(size=12)) +
                                                           labs(fill = glue("Mean Absolute T-Statistic\nfor {SPI_nickname}")) +
                                                           guides(fill = guide_colorbar(title.position = "top", 
                                                                                        nrow = 1,
                                                                                        barwidth = 12, 
                                                                                        barheight = 0.75,
                                                                                        title.hjust = 0.5,
                                                                                        label.position = "bottom"))))
      
      # Append to list
      dataset_ggseg <- dataset_ggseg + theme(legend.position = "none")
      ggseg_plot_list <- list.append(ggseg_plot_list, dataset_ggseg)
      ggseg_plot_list <- list.append(ggseg_plot_list, plot_spacer())
      legend_list <- list.append(legend_list, ctx_legend)
    }
  }
  
  output <- list(brain_plots = ggseg_plot_list,
                 legend_plots = legend_list)
  return(output)
}
  

# Pearson correlation in the brain -- mean abs t-statistic by region across all connections
pearson_regional_data <- pairwise_t_stats_by_region_from %>%
  filter(SPI == "cov_EmpiricalCovariance") %>%
  left_join(., UCLA_CNP_brain_region_info)  %>%
  dplyr::select(-Index) %>%
  left_join(., ABIDE_ASD_brain_region_info)
correlation_plots <- plot_data_per_group(SPI_nickname = "Pearson R",
                    input_data = pearson_regional_data,
                    group_colors = group_colors)

wrap_plots(correlation_plots[[1]], 
           ncol=2, 
           byrow=T)
ggsave(glue("{plot_path}/Region_wise_avg_t_stat_pearson_corrs.png"),
       width=4, height=7, units="in", dpi=300)
wrap_plots(correlation_plots[[2]], 
           nrow=4, 
           byrow=T)
ggsave(glue("{plot_path}/Region_wise_avg_t_stat_pearson_corrs_legend.png"),
       width=3, height=4, units="in", dpi=300)

# Gaussian directed information: FROM each region
min_max_df <- pairwise_t_stats_by_region_from %>%
  plyr::rbind.fill(pairwise_t_stats_by_region_to) %>%
  filter(SPI == "di_gaussian") %>%
  group_by(Study, Comparison_Group) %>%
  summarise(min_fill = min(mean_T_magnitude),
            max_fill = max(mean_T_magnitude))

# Get min per group
DI_from_regional_data <- pairwise_t_stats_by_region_from %>%
  filter(SPI == "di_gaussian") %>%
  left_join(., UCLA_CNP_brain_region_info)  %>%
  dplyr::select(-Index) %>%
  left_join(., ABIDE_ASD_brain_region_info)
DI_from_plots <- plot_data_per_group(SPI_nickname = "Directed Information",
                                         input_data = DI_from_regional_data,
                                         group_colors = group_colors,
                                     min_max_df = min_max_df)
wrap_plots(DI_from_plots[[1]], 
           ncol=2, 
           byrow=T)
ggsave(glue("{plot_path}/Region_wise_avg_t_stat_gaussian_DI_from.png"),
       width=4, height=7, units="in", dpi=300)
wrap_plots(DI_from_plots[[2]], 
           nrow=4, 
           byrow=T)
ggsave(glue("{plot_path}/Region_wise_avg_t_stat_gaussian_DI_from_legend.png"),
       width=3, height=4, units="in", dpi=300)

# Gaussian directed information: TO each region
DI_to_regional_data <- pairwise_t_stats_by_region_to %>%
  filter(SPI == "di_gaussian") %>%
  left_join(., UCLA_CNP_brain_region_info)  %>%
  dplyr::select(-Index) %>%
  left_join(., ABIDE_ASD_brain_region_info)
DI_to_plots <- plot_data_per_group(SPI_nickname = "Directed Information",
                                         input_data = DI_to_regional_data,
                                         group_colors = group_colors,
                                   min_max_df = min_max_df)
wrap_plots(DI_to_plots[[1]], 
           ncol=2, 
           byrow=T)
ggsave(glue("{plot_path}/Region_wise_avg_t_stat_gaussian_DI_to.png"),
       width=4, height=7, units="in", dpi=300)
wrap_plots(DI_to_plots[[2]], 
           nrow=4, 
           byrow=T)
ggsave(glue("{plot_path}/Region_wise_avg_t_stat_gaussian_DI_to_legend.png"),
       width=3, height=4, units="in", dpi=300)



# Spaghetti plot to vs from for gaussian DI per brain region by group
pairwise_t_stats_by_region_from %>%
  plyr::rbind.fill(pairwise_t_stats_by_region_to) %>%
  filter(SPI == "di_gaussian") %>%
  mutate(Comparison_Group = factor(Comparison_Group, levels = c("SCZ", "BPD", "ADHD", "ASD")))%>%
  ggplot(data=., mapping=aes(x=Direction, y=mean_T_magnitude, group=Brain_Region, color=Comparison_Group)) +
  geom_line(linewidth=0.4, alpha=0.5) +
  facet_wrap(. ~ Comparison_Group, scales="free_y", nrow=1) +
  ylab("Mean T-Statistic for DI\nby Brain Region") +
  xlab("Functional Connectivity Direction")  +
  scale_color_manual(values = c("SCZ"="#573DC7", 
                                "BPD"="#D5492A", 
                                "ADHD"="#0F9EA9",
                                "ASD"="#C47B2F"))  +
  theme(legend.position = "none")
ggsave(glue("{plot_path}/DI_Gaussian_T_stats_to_vs_from.png"),
       width=10, height=2.5, units="in", dpi=300)


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