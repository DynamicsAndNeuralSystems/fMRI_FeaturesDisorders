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
this_plot <- plot_feature_in_brain(fill_color_gradient=brain_colors[12], 
                      region_label="all")
this_plot
ggsave(glue("{plot_path}/Example_univariate_feature_brain.svg"),
       width=3, height=2, units="in", dpi=300)

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

# Wang's periodicity infographic

# Isolate brain regions of interest
if (!(file.exists(glue("{data_path}/BOLD_data_for_wangs_PD.csv")))) {
  UCLA_CNP_metadata <- pyarrow_feather$read_feather("~/data/UCLA_CNP/study_metadata/UCLA_CNP_sample_metadata.feather") %>%
    mutate(Study="UCLA_CNP")
  
  UCLA_CNP_catch24 <- pyarrow_feather$read_feather("~/data/UCLA_CNP/processed_data/UCLA_CNP_AROMA_2P_GMR_catch24_filtered.feather")  %>%
    left_join(., UCLA_CNP_sample_metadata)
  
  Raw_TS_data_UCLA_CNP <- pyarrow_feather$read_feather("~/data/UCLA_CNP/raw_data/UCLA_CNP_AROMA_2P_GMR_fMRI_TS.feather")
  
  # 
  compute_ACF_for_region <- function(split_df) {
    region_ACF <- acf(split_df$values, 
                      lag.max = (1/3)*length(unique(split_df$timepoint)),
                      type = "correlation",
                      plot = FALSE, na.action = na.fail, demean = TRUE)
    region_ACF_df <- data.frame(Brain_Region = unique(split_df$Brain_Region),
                                Sample_ID = unique(split_df$Sample_ID),
                                ACF = region_ACF$acf,
                                Lag = region_ACF$lag - 1) %>%
      filter(Lag > 0)
    return(region_ACF_df)
  }
  TS_data_broken_down <- Raw_TS_data_UCLA_CNP %>%
    group_by(Sample_ID, Brain_Region) %>%
    group_split()
  
  All_samples_regions_ACF <- do.call(plyr::rbind.fill, TS_data_broken_down %>%
                                       purrr::map(~ compute_ACF_for_region(split_df=.x)))
  
  # Filter to only lags +/- periodicity value and only ACF > 0.1
  All_samples_regions_ACF %>%
    left_join(., UCLA_CNP_catch24 %>% filter(names=="PD_PeriodicityWang_th0_01")) %>%
    group_by(Sample_ID, Brain_Region) %>%
    filter(between(Lag, values-1, values+1)) %>%
    filter(n() == 3) %>%
    mutate(diff_before = ACF[Lag==values] - ACF[Lag < values],
           diff_after = ACF[Lag==values] - ACF[Lag > values]) %>%
    ungroup() %>%
    filter(ACF > 0.15, Lag==values, Lag > 2) %>%
    filter(diff_before > 0.1, diff_after > 0.1) %>%
    group_by(Sample_ID) %>%
    filter(Lag %in% c(min(Lag), max(Lag))) %>%
    mutate(Lag_Diff = max(Lag) - min(Lag)) %>%
    arrange(desc(Lag_Diff))
  
  BOLD_data_for_wangs_PD <- Raw_TS_data_UCLA_CNP  %>%
    filter(Sample_ID == "sub-10440",
           Brain_Region %in% c("ctx-lh-caudalmiddlefrontal",
                               "ctx-lh-parstriangularis"))
  
  write.table(BOLD_data_for_wangs_PD,
              glue("{data_path}/BOLD_data_for_wangs_PD.csv"),
              sep=",", row.names = F)
} else {
  BOLD_data_for_wangs_PD <- read.csv(glue("{data_path}/BOLD_data_for_wangs_PD.csv"))
}

# Read in cubic spline detrended time series from Matlab
BOLD_data_for_wangs_PD_DT <- read.csv(glue("{data_path}/brain_region_TS_detrended.csv"))


# Low periodicity TS
BOLD_data_for_wangs_PD %>%
  filter(Brain_Region == "ctx-lh-caudalmiddlefrontal") %>%
  ggplot(data=., mapping=aes(x=timepoint, y=values)) +
  geom_line(color = "blue") +
  xlab("Time (s)") +
  ylab("Value")
ggsave(glue("{plot_path}/Periodicity_Low_TS.svg"), width=3, height=1.25, units="in", dpi=300)

#
low_periodicity_ACF <- acf(BOLD_data_for_wangs_PD_DT %>%
                             filter(Periodicity_Type=="Low") %>% 
                             pull(Detrended_Value), 
                           lag.max = (1/3)*length(unique(BOLD_data_for_wangs_PD_DT$Timepoint)),
                           type = "correlation",
                           plot = FALSE, na.action = na.fail, demean = TRUE)
data.frame(Lag=low_periodicity_ACF$lag,
           ACF=low_periodicity_ACF$acf) %>%
  mutate(Lag = Lag-1) %>%
  filter(Lag>0) %>%
  ggplot(data=., mapping = aes(x=Lag,
                               y=ACF)) + 
  geom_line() +
  geom_vline(xintercept = 3, color="blue")
ggsave(glue("{plot_path}/Periodicity_Low_ACF.svg"), width=3, height=1.25, units="in", dpi=300)

# High periodicity TS
BOLD_data_for_wangs_PD %>%
  filter(Brain_Region == "ctx-lh-parstriangularis") %>%
  ggplot(data=., mapping=aes(x=timepoint, y=values)) +
  geom_line(color = "red") +
  xlab("Time (s)") +
  ylab("Value")
ggsave(glue("{plot_path}/Periodicity_High_TS.svg"), width=3, height=1.25, units="in", dpi=300)

# Compute autocorrelation for each TS
high_periodicity_ACF <- acf(BOLD_data_for_wangs_PD_DT %>%
                             filter(Periodicity_Type=="High") %>% 
                             pull(Detrended_Value), 
                           lag.max = (1/3)*length(unique(BOLD_data_for_wangs_PD_DT$Timepoint)),
                           type = "correlation",
                           plot = FALSE, na.action = na.fail, demean = TRUE)
data.frame(Lag = high_periodicity_ACF$lag,
           ACF = high_periodicity_ACF$ac) %>%
  mutate(Lag = Lag-1) %>%
  filter(Lag > 0) %>%
  ggplot(data=., mapping = aes(x=Lag,
                               y=ACF)) + 
  geom_line() +
  geom_vline(xintercept = 23, color="red")
ggsave(glue("{plot_path}/Periodicity_High_ACF.svg"), width=3, height=1.25, units="in", dpi=300)


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
ggsave(glue("{plot_path}/Pairwise_feature_brains.svg"),
       width=9, height=5, units="in", dpi=300, bg="white")


################################################################################
# Pearson correlation infographic
# Original time series
sr <- 100 # sampling rate
ts <- 1.0/sr # sampling interval
t <- seq(0, 1, ts)
freq <- 1
x <- 2*sin(2*pi*freq*t)
freq <- 4
x <- x + sin(2*pi*freq*t)
freq <- 7
x <- x + 0.75*sin(2*pi*freq*t)

# Two time series with high magnitude correlation
set.seed(127)
x_corr1 <-  x + 0.9*rnorm(101)
x_corr2 <- x + 0.9*rnorm(101)

# Line plot for the two time series
data.frame(TS1 = x_corr1,
           TS2 = x_corr2,
           timepoint = 1:101) %>%
  pivot_longer(cols=c(TS1,TS2)) %>%
  ggplot(data=., mapping=aes(x=timepoint, y=value, group=name, color=name)) +
  geom_line() +
  xlab("Time (s)") +
  ylab("Value") +
  theme(legend.position="none") +
  scale_color_manual(values=c("#f88c4d", "#4db9f8"))
ggsave(glue("{plot_path}/TS_Pair_HighCorr.svg"), width=3, height=1.25, units="in", dpi=300)

data.frame(TS1 = x_corr1,
           TS2 = x_corr2,
           timepoint = 1:101) %>%
  ggplot(data=., mapping=aes(x=TS1, y=TS2)) +
  geom_point(color="black", size=0.1) +
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks = element_blank()) +
  coord_equal()
ggsave(glue("{plot_path}/TS_Pair_HighCorr_Points.svg"), width=1, height=1, units="in", dpi=300)

cor(x_corr1,x_corr2, method="pearson")


# Two time series with low magnitude correlation
set.seed(127)
x_corr1 <-  x + 0.9*rnorm(101)
x_corr2 <- rnorm(101, mean=0, sd=2)

# Line plot for the two time series
data.frame(TS1 = x_corr1,
           TS2 = x_corr2,
           timepoint = 1:101) %>%
  pivot_longer(cols=c(TS1,TS2)) %>%
  ggplot(data=., mapping=aes(x=timepoint, y=value, group=name, color=name)) +
  geom_line() +
  xlab("Time (s)") +
  ylab("Value") +
  theme(legend.position="none") +
  scale_color_manual(values=c("#f88c4d", "#4db9f8"))
ggsave(glue("{plot_path}/TS_Pair_LowCorr.svg"), width=3, height=1.25, units="in", dpi=300)

data.frame(TS1 = x_corr1,
           TS2 = x_corr2,
           timepoint = 1:101) %>%
  ggplot(data=., mapping=aes(x=TS1, y=TS2)) +
  geom_point(color="black", size=0.1) +
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks = element_blank()) +
  scale_x_continuous(limits=c(min(x_corr1), max(x_corr1))) +
  scale_y_continuous(limits=c(min(x_corr1), max(x_corr1))) +
  coord_equal()
ggsave(glue("{plot_path}/TS_Pair_LowCorr_Points.svg"), width=1, height=1, units="in", dpi=300)

cor(x_corr1,x_corr2, method="pearson")
