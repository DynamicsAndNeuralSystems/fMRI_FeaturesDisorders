################################################################################
# Define study/data paths
################################################################################

github_dir <- "~/github/fMRI_FeaturesDisorders/"
plot_path <- paste0(github_dir, "plots/Manuscript_Draft/methods_overview/")
TAF::mkdir(plot_path)

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
ggsave(glue("{plot_path}/example_brain.png"),
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
ggsave(glue("{plot_path}/example_brain_univariate.png"),
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
ggsave(glue("{plot_path}/example_brain_pairwise.png"),
       width=3,height=2,units="in", dpi=300)

# Sample feature vector
data.frame(x=sample(1:24),
           y=as.character(1:24)) %>%
  ggplot(data=., mapping=aes(x=x, y=0, fill=y)) +
  geom_tile(color="black") +
  scale_fill_viridis_d() +
  theme_void() +
  theme(legend.position = "none")
ggsave(glue("{plot_path}/example_feature_vector.png"),
       width=3,height=0.2,units="in", dpi=300)

# Brains with 24 example features
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

brain_colors <- viridis::viridis(24)

# Plot the full brain
plots <- sample(1:24) %>%
  purrr::map(~ plot_feature_in_brain(fill_color_gradient=brain_colors[.x], 
                                     region_label="all"))
wrap_plots(plots, nrow=4)
ggsave(glue("{plot_path}/Univariate_feature_brains.png"),
       width=9, height=4.5, units="in", dpi=300, bg="white")


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
ggsave(glue("{plot_path}/example_pyspi_feature_vector.png"),
       width=3,height=0.2,units="in", dpi=300)

# Example network graph
set.seed(127)

# Data for all pairwise connection graphs
# Edges are defined as cortical lobe --> specific ROI connection
edges <- read.csv("region_node_hierarchy.csv") %>% distinct()
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
                    size=2) +
    scale_edge_color_gradientn(colors=c(alpha(edge_color,0.3), edge_color)) +
    theme_void() + 
    theme(plot.title=element_text(size=14, face="bold", hjust=0.5),
          legend.position="none")
  
  return(p)
}

plots <- sample(1:14) %>%
  purrr::map(~ plot_network_data(pyspi_colors[.x]))
wrap_plots(plots, nrow=3)
ggsave(glue("{plot_path}/Pairwise_feature_brains.png"),
       width=9, height=5, units="in", dpi=300, bg="white")