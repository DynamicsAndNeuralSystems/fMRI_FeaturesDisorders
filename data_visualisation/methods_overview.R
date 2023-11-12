################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/Library/CloudStorage/OneDrive-TheUniversityofSydney(Students)/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/methods_overview/")
TAF::mkdir(plot_path)
data_path <- "~/data/TS_feature_manuscript"

library(tidyverse)
library(glue)
library(icesTAF)
library(cowplot)
library(ggseg)
library(patchwork)
library(ggseg)
library(LaCroixColoR)
library(tidyverse)
library(ggraph)
library(igraph)
library(dendextend)
theme_set(theme_cowplot())

python_to_use <- "/Users/abry4213/anaconda3/envs/pyspi/bin/python3"
reticulate::use_python(python_to_use)

library(reticulate)

# Import pyarrow.feather as pyarrow_feather
pyarrow_feather <- import("pyarrow.feather")

################################################################################
# fMRI --> MTS overview
################################################################################

ggplot() +
  geom_brain(atlas = dk, hemi="left", side="lateral", color="gray30") +
  theme_void() + 
  theme(legend.position="none")
ggsave(glue("{plot_path}/ggseg_full_color.svg"),
       width=3,height=2,units="in", dpi=300)

# Plot a few example time-series
example_MTS <- as.data.frame(matrix(data=rnorm(500), nrow=50, ncol=10)) %>%
  mutate(Timepoint = 1:50) %>%
  pivot_longer(cols=c(-Timepoint), names_to = "Brain_Region", values_to="BOLD") %>%
  mutate(Brain_Region = as.factor(as.numeric(gsub("V", "", Brain_Region))))

# Lines
example_MTS %>%
  ggplot(mapping=aes(x=Timepoint, y=BOLD, color=Brain_Region)) +
  geom_line() + 
  facet_grid(Brain_Region ~ .) +
  theme_void() +
  theme(legend.position = "none",
        strip.text = element_blank(),
        panel.spacing = unit(-1, "lines"))
ggsave(glue("{plot_path}/ggseg_example_TS_lines.svg"),
       width=3,height=2,units="in", dpi=300)

# MTS heatmap
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

MTS_heatmap_list <- list()
MTS_colors <- gg_color_hue(10)
for (i in 1:length(unique(example_MTS$Brain_Region))) {
  roi = unique(example_MTS$Brain_Region)[i]
  roi_color = MTS_colors[i]
  roi_p <- example_MTS %>%
    filter(Brain_Region == roi) %>%
    ggplot(mapping=aes(x=Timepoint, fill=BOLD, y=0)) +
    geom_tile() + 
    theme_void() +
    scale_fill_gradient(low=alpha(roi_color, 0.2), high=roi_color) +
    theme(legend.position = "none")
  MTS_heatmap_list[[roi]] <- roi_p
}
wrap_plots(MTS_heatmap_list, ncol=1)
ggsave(glue("{plot_path}/ggseg_example_TS_heatmap.svg"),
       width=3,height=2,units="in", dpi=300)
  
################################################################################
# Univariate methods
################################################################################

# Example brain with regions highlighted
dk %>%
  as_tibble() %>%
  mutate(fillval = case_when(label == "lh_caudalmiddlefrontal" ~ "1",
                             label == "lh_insula" ~ "2",
                             label == "lh_lateraloccipital" ~ "3", 
                             T ~ NA_character_)) %>%
  ggseg(atlas = "dk", mapping = aes(fill = fillval),
        hemisphere="left",
        view = "lateral",
        position = "stacked", colour = "gray50") +
  scale_fill_manual(values=c("#ED702D", "#4E70BF", "#B12718"),
                    na.value="white") +
  theme_void() +
  theme(plot.title = element_blank(),
        legend.position = "none")
ggsave(glue("{plot_path}/example_brain.svg"),
       width=3,height=2,units="in", dpi=300)

# Example brain, univariate
dk %>%
  as_tibble() %>%
  mutate(fillval = case_when(label == "lh_caudalmiddlefrontal" ~ "1",
                             T ~ NA_character_)) %>%
  ggseg(atlas = "dk", mapping = aes(fill = fillval),
        hemisphere="left",
        view = "lateral",
        position = "stacked", colour = "gray50") +
  scale_fill_manual(values=c("#ED702D"),
                    na.value="white") +
  theme_void() +
  theme(plot.title = element_blank(),
        legend.position = "none")
ggsave(glue("{plot_path}/example_brain_univariate.svg"),
       width=3,height=2,units="in", dpi=300)

# Example brain, pairwise
dk %>%
  as_tibble() %>%
  mutate(fillval = case_when(label == "lh_insula" ~ "1",
                             label == "lh_lateraloccipital" ~ "2", 
                             T ~ NA_character_)) %>%
  ggseg(atlas = "dk", mapping = aes(fill = fillval),
        hemisphere="left",
        view = "lateral",
        position = "stacked", colour = "gray50") +
  scale_fill_manual(values=c("#4E70BF", "#B12718"),
                    na.value="white") +
  theme_void() +
  theme(plot.title = element_blank(),
        legend.position = "none")
ggsave(glue("{plot_path}/example_brain_pairwise.svg"),
       width=3,height=2,units="in", dpi=300)

# Sample feature vector
data.frame(x=sample(1:24),
           y=as.character(1:24)) %>%
  ggplot(data=., mapping=aes(x=x, y=0, fill=y)) +
  geom_tile(color="black") +
  scale_fill_viridis_d() +
  theme_void() +
  theme(legend.position = "none")
ggsave(glue("{plot_path}/example_feature_vector.svg"),
       width=3,height=0.2,units="in", dpi=300)


brain_colors <- viridis::viridis(25)

# Sample region vector
data.frame(x=sample(1:35),
           y=as.character(1:35)) %>%
  ggplot(data=., mapping=aes(x=y, y=0, fill=x)) +
  geom_tile(color="black") +
  scale_fill_gradientn(colors=c(alpha(brain_colors[12], 0.3), 
                                brain_colors[12]), 
                       na.value=NA) +
  theme_void() +
  theme(legend.position = "none")
ggsave(glue("{plot_path}/example_region_vector.svg"),
       width=3,height=0.2,units="in", dpi=300)

# Brains with 25 example features
plot_feature_in_brain <- function(fill_color_gradient, region_label="all") {
  if (region_label=="all") {
    p <- dk %>%
      as_tibble() %>%
      mutate(region_values = runif(nrow(.))) %>%
      ggseg(atlas = "dk", mapping = aes(fill = region_values),
            hemisphere="left",
            view = "lateral",
            position = "stacked", colour = "black") +
      scale_fill_gradientn(colors=c(alpha(fill_color_gradient, 0.3), 
                                    fill_color_gradient), 
                           na.value=NA)
  } else {
    p <- dk %>%
      as_tibble() %>%
      mutate(region_values = ifelse(label==region_label, "1", NA_character_)) %>%
      ggseg(atlas = "dk", mapping = aes(fill = region_values),
            hemisphere="left",
            view = "lateral",
            position = "stacked", colour = "gray40") +
      scale_fill_manual(values=c(fill_color_gradient),
                        na.value="white")
  }

  p <- p  +
    theme_void() +
    theme(plot.title = element_blank(),
          legend.position = "none") 
}


# Plot the full brain
plots <- sample(1:25) %>%
  purrr::map(~ plot_feature_in_brain(fill_color_gradient=brain_colors[.x], 
                                     region_label="all"))
wrap_plots(plots, nrow=5)
ggsave(glue("{plot_path}/Univariate_feature_brains.svg"),
       width=6, height=5, units="in", dpi=300, bg="white")

# Plot just one example feature
plot_feature_in_brain(fill_color_gradient="#548AC5", 
                      region_label="all")
ggsave(glue("{plot_path}/Example_univariate_feature_brain.svg"),
       width=3, height=2, units="in", dpi=300)
data.frame(x=sample(1:10),
           y=as.character(1:10)) %>%
  ggplot(data=., mapping=aes(y=y, x=0, fill=x)) +
  geom_tile(color="black") +
  scale_fill_gradientn(colors=c(alpha("#548AC5", 0.3), 
                                "#548AC5"), 
                       na.value=NA) +
  theme_void() +
  theme(legend.position = "none")
ggsave(glue("{plot_path}/Example_univariate_feature_brain_vector.svg"),
       width=1, height=2, units="in", dpi=300)

# Plot feature X region matrix
feature_vec_plot_list <- list()
for (i in sample(1:25)) {
  p <- data.frame(x=sample(1:10),
             y=as.character(1:10)) %>%
    ggplot(data=., mapping=aes(x=y, y=0, fill=x)) +
    geom_raster(color="black") +
    scale_fill_gradientn(colors=c(alpha(brain_colors[i], 0.3), 
                                  brain_colors[i]), 
                         na.value=NA) +
    theme_void() +
    theme(legend.position = "none")
  feature_vec_plot_list <- list.append(feature_vec_plot_list, p)
}

wrap_plots(feature_vec_plot_list, ncol=1)
ggsave(glue("{plot_path}/example_region_feature_vector.svg"),
       width=2,height=5,units="in", dpi=300)

################################################################################

################################################################################
# Pairwise methods
################################################################################

pyspi_colors <- unname(lacroix_palette("PassionFruit",14))

# Sample feature vector
data.frame(x=sample(1:14),
           y=as.character(1:14)) %>%
  ggplot(data=., mapping=aes(x=x, y=0, fill=y)) +
  geom_tile(color="black") +
  scale_fill_manual(values=pyspi_colors) +
  theme_void() +
  theme(legend.position = "none")
ggsave(glue("{plot_path}/example_pyspi_feature_vector.svg"),
       width=3,height=0.2,units="in", dpi=300)

# Example network graph
set.seed(127)

# Data for all pairwise connection graphs
# Edges are defined as cortical lobe --> specific ROI connection
edges <- read.csv(glue("{github_dir}/data_visualisation/region_node_hierarchy.csv")) %>% distinct()
# ROIs don't include the origin --> cortical lobe connection
rois <- edges %>% filter(!(to %in% c("Cingulate", "Frontal", "Insula",
                                     "Occipital", "Parietal", "Temporal")))
# Create a dataframe of vertices, one line per object in the ROI cortical lobe hierarchy
vertices = data.frame(name = unique(c(as.character(edges$from), as.character(edges$to))))
vertices$group <- edges$from[match(vertices$name, edges$to)]
# Create an igraph object
mygraph <- graph_from_data_frame(edges, vertices=vertices)

# Function to plot specific data in a given iteration
plot_network_data <- function(edge_color) {
  
  # connect = dataframe of pairwise correlations between cortical ROIs
  connect <- data.frame("from" = sample(rois$to, 500, replace=T),
                        "to" = sample(rois$to,500, replace=T),
                        value = runif(500, 0, 1)) %>%
    arrange(from, to) %>%
    sample_n(100)
  
  # Collect edges where from = selected ROI and to = ROIs connected to the selected ROI
  from <- match(connect$from, vertices$name)
  to <- match(connect$to, vertices$name)
  
  # mygraph = igraph object linking each cortical ROI
  # convert to a circular dendrogram-shaped ggraph object
  p <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
    theme_void() +  geom_conn_bundle(data = get_con(from = from, to = to, 
                                                    value=connect$value), 
                                     tension=runif(1, 0.55, 0.85), 
                                     width=1,
                                     aes(color=value,
                                         alpha=value))  +
    labs(edge_width="Pearson\nCorrelation") + 
    geom_node_point(aes(filter = leaf, 
                        x = x*1.05, y=y*1.05),   
                    size=3) +
    scale_edge_color_gradientn(colors=c(alpha(edge_color,0.3), edge_color)) +
    theme_void() + 
    theme(plot.title=element_text(size=14, face="bold", hjust=0.5),
          legend.position="none")
  
  return(p)
}

plots <- sample(1:14) %>%
  purrr::map(~ plot_network_data(pyspi_colors[.x]))
wrap_plots(plots, nrow=3)
ggsave(glue("{plot_path}/Pairwise_feature_brains.svg"),
       width=9, height=5, units="in", dpi=300, bg="white")


# Plot example purple network
plot_network_data(edge_color="#BF91D0")
ggsave(glue("{plot_path}/Example_network.svg"),
       width=3, height=3, units="in", dpi=300, bg="white")
