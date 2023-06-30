import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import nxviz as nv
from nxviz import lines, utils, layouts, nodes, plots, edges, encodings
from random import choice
from nxviz import encodings as aes
from typing import Callable, Dict, Hashable, Tuple, Optional
from matplotlib.colors import ListedColormap, Normalize, BoundaryNorm
from functools import partial
from matplotlib.cm import get_cmap
import pylab as pl
from sklearn import preprocessing
import matplotlib.ticker as ticker

###############################################################################
# Helper functions
# Function to color nodes by group
def group_colormap(data: pd.Series):
    cmap = {"Cingulate": "#A6CEE3", 
            "Frontal": "#2078B4", 
            "Insula": "#B2DF8A",
            "Occipital": "#31A02D",
            "Parietal": "#FB9A98",
            "Subcortex": "#E21A1B",
            "Temporal": "#FDBF6E"}
    return data.apply(lambda x: cmap.get(x))
    
# Custom edge coloring functions
def continuous_color_func(val, cmap, data: pd.Series, vmin, vmax):
    """Return RGBA of a value.

    ## Parameters

    - `val`: Value to convert to RGBA
    - `cmap`: A Matplotlib cmap
    - `data`: Pandas series.
    """
    norm = Normalize(vmin=vmin, vmax=vmax)
    return cmap(norm(val))

def color_func(data: pd.Series, cmap, vmin, vmax) -> Callable:
    """Return a color function that takes in a value and returns an RGB(A) tuple.

    This will do the mapping to the continuous and discrete color functions.
    """
    func = continuous_color_func
    return partial(func, cmap=cmap, data=data, vmin=vmin, vmax=vmax)

def data_color(data: pd.Series, ref_data: pd.Series, cmap, vmin, vmax) -> pd.Series:
    """Return iterable of colors for a given data.

    `cfunc` gives users the ability to customize the color mapping of a node.
    The only thing that we expect is that it takes in a value
    and returns a matplotlib-compatible RGB(A) tuple or hexadecimal value.

    The function takes in `ref_data`
    which is used to determine important colormap values (such as boundaries).
    That colormap is then applied to the actual `data`.

    ## Parameters

    - `data`: The data on which to map colors.
    - `ref_data`: The data on which the colormap is constructed.
    """
    cfunc = color_func(ref_data, cmap, vmin, vmax)
    return data.apply(cfunc)

def edge_colors(
    et: pd.DataFrame,
    nt: pd.DataFrame,
    cmap,
    vmin,
    vmax,
    color_by: Hashable,
    node_color_by: Hashable,
):
    """Default edge line color function."""
    if color_by in ("source_node_color", "target_node_color"):
        edge_select_by = color_by.split("_")[0]
        return data_color(
            et[edge_select_by].apply(nt[node_color_by].get),
            nt[node_color_by], cmap=cmap, vmin=vmin, vmax=vmax)
    elif color_by:
        return data_color(et[color_by], et[color_by], cmap=cmap, vmin=vmin, vmax=vmax)
    return pd.Series(["black"] * len(et), name="color_by")


###############################################################################

# Load general lobe-region hierarchy df
node_to_from_structure =  pd.read_csv("/Users/abry4213/data/TS_feature_manuscript/node_to_from_structure.csv")

# Define a formatting function to add a minus sign to the hypo-connectivity labels
def format_func(value, tick_number):
    return f'-{value:.1f}'

def plot_individual_chord_graph(G, my_cmap, vmin, vmax):
    # Customize node styling
    ax = plt.gca()

    nt = utils.node_table(G)
    pos = layouts.circos(nt, group_by="lobe")

    # Define nodes and corresponding lobe colors
    node_color = group_colormap(nt["lobe"])
    alpha = nodes.transparency(nt, alpha_by=None)
    size = nodes.node_size(nt, size_by=None)
    patches = nodes.node_glyphs(
        nt, pos, node_color=node_color, alpha=alpha, size=size)
    for patch in patches:
        ax.add_patch(patch)
        
    # Customize edge styling
    et = utils.edge_table(G)
    edge_color = edge_colors(et, nt=None, color_by="weight", node_color_by=None, cmap=my_cmap, vmin=vmin, vmax=vmax)
    lw = pd.Series([2]*et.shape[0])
    edge_alpha = pd.Series(preprocessing.minmax_scale(np.absolute(et["weight"]), feature_range=(0.4, 1), axis=0, copy=True))
    patches = lines.circos(
        et, pos, edge_color=edge_color, alpha=edge_alpha, lw=lw, aes_kw={"fc": "none"}
    )
    for patch in patches:
        ax.add_patch(patch)
        

    plots.rescale(G)
    plots.aspect_equal()
    plots.despine()

def plot_hyper_hypo_chord_per_group(study, comparison_group, data_path,
                              plot_path,
                              node_to_from_structure,
                              hyper_threshold = 0.6, 
                              hypo_threshold = -0.6, 
                              hyper_cmap = "YlOrRd", 
                              hypo_cmap = "YlGnBu"):
    
    # Load all connections for group
    group_FC_connections = pd.read_csv(f"{data_path}/{study}_{comparison_group}_all_FC_pairs.csv") 
    group_node_lookup = node_to_from_structure.query("Study == @study")
    
    # Convert connections to an array
    group_FC_arr = np.array(group_FC_connections)
    
    # Create one directed graph for hyper and hypo connectivity
    G_hyper = nx.DiGraph(name = "Hyperconnectivity")
    G_hypo = nx.DiGraph(name = "Hypo")
    
    # Add all nodes
    G_hyper.add_nodes_from(group_node_lookup.to.to_list())
    G_hypo.add_nodes_from(group_node_lookup.to.to_list())
    
    # Iterate over the nodes and add weights above or below the specified thresholds only
    for i in range(len(group_FC_arr)):
        connection = group_FC_arr[i]
        region_from = connection[0] # protein a node
        region_to = connection[1] # protein b node
        mean_beta = connection[2] # score as weighted edge
        
        # Add only high edges to G_hyper
        if mean_beta > hyper_threshold:
            G_hyper.add_weighted_edges_from([(region_from,region_to,mean_beta)])
        # Add only low edges to G_hypo
        if mean_beta < hypo_threshold:
            G_hypo.add_weighted_edges_from([(region_from,region_to,-1*mean_beta)])
            
    # Add lobe to each node
    for n, d in G_hyper.nodes(data=True):
        brain_region = n
        lobe = node_to_from_structure.query("to == @brain_region")["from"].to_string(index=False)
        G_hyper.nodes[n]["lobe"] = lobe
    for n, d in G_hypo.nodes(data=True):
        brain_region = n
        lobe = node_to_from_structure.query("to == @brain_region")["from"].to_string(index=False)
        G_hypo.nodes[n]["lobe"] = lobe
        
    # Figure out edge color ranges to match
    G_hyper_et = utils.edge_table(G_hyper)
    G_hypo_et = utils.edge_table(G_hypo)
    
    vmin = min(G_hypo_et.weight.tolist() + G_hyper_et.weight.tolist())
    vmax = max(G_hypo_et.weight.tolist() + G_hyper_et.weight.tolist())
    a = np.array([[vmin,vmax]])
        
    ###########################################################################
    # Hyper plot
    fig_hyper, ax_hyper = plt.subplots()
    plot_individual_chord_graph(G_hyper, get_cmap(hyper_cmap), vmin = vmin, vmax = vmax)
    fig_hyper.tight_layout()
    fig_hyper.set_size_inches(3,3)
    fig_hyper.savefig(f"{plot_path}/{study}_{comparison_group}_hyper_FC.svg", dpi=300,
               bbox_inches='tight', pad_inches=0.0, transparent=True)

    # Hyper legend
    pl.figure(figsize=(0.5, 5))
    img = pl.imshow(a, cmap=hyper_cmap)
    pl.gca().set_visible(False)
    cax = pl.axes([0.1, 0.2, 0.8, 0.6])
    hyper_cbar = pl.colorbar(cax=cax)
    # Increase the font size of the colorbar labels
    hyper_cbar.ax.yaxis.set_tick_params(labelsize=20)
    pl.savefig(f"{plot_path}/{study}_{comparison_group}_hyper_FC_legend.svg", dpi=300,
               bbox_inches='tight', pad_inches=0.0, transparent=True)

    # Hypo plot
    fig_hypo, ax_hypo = plt.subplots()
    plot_individual_chord_graph(G_hypo, get_cmap(hypo_cmap), vmin = vmin, vmax = vmax)
    fig_hypo.tight_layout()
    fig_hypo.set_size_inches(3,3)
    fig_hypo.savefig(f"{plot_path}/{study}_{comparison_group}_hypo_FC.svg", dpi=300,
               bbox_inches='tight', pad_inches=0.0, transparent=True)
    
    # Hypo legend
    pl.figure(figsize=(0.5, 5))
    img = pl.imshow(a, cmap=hypo_cmap)
    pl.gca().set_visible(False)
    cax = pl.axes([0.1, 0.2, 0.8, 0.6])
    hypo_cbar = pl.colorbar(cax=cax)
    hypo_cbar.ax.yaxis.set_major_formatter(ticker.FuncFormatter(format_func))
    
    # Increase the font size of the colorbar labels
    hypo_cbar.ax.yaxis.set_tick_params(labelsize=20)
    pl.savefig(f"{plot_path}/{study}_{comparison_group}_hypo_FC_legend.svg", dpi=300,
               bbox_inches='tight', pad_inches=0.0, transparent=True)


###############################################################################
# Call for UCLA CNP SCZ
plot_hyper_hypo_chord_per_group(study = "UCLA_CNP", 
                          comparison_group = "SCZ", 
                          data_path = "/Users/abry4213/data/TS_feature_manuscript/",
                          plot_path = "/Users/abry4213/github/fMRI_FeaturesDisorders/plots/Manuscript_Draft/pairwise_results/",
                          node_to_from_structure = node_to_from_structure,
                            hyper_threshold = 0.6, 
                            hypo_threshold = -0.6, 
                            hyper_cmap = "Reds", 
                            hypo_cmap = "Blues")

# UCLA CNP BPD
plot_hyper_hypo_chord_per_group(study = "UCLA_CNP", 
                          comparison_group = "BPD", 
                          data_path = "/Users/abry4213/data/TS_feature_manuscript/",
                          plot_path = "/Users/abry4213/github/fMRI_FeaturesDisorders/plots/Manuscript_Draft/pairwise_results/",
                          node_to_from_structure = node_to_from_structure,
                            hyper_threshold = 0.6, 
                            hypo_threshold = -0.6, 
                            hyper_cmap = "Reds", 
                            hypo_cmap = "Blues")

# UCLA CNP ADHD
plot_hyper_hypo_chord_per_group(study = "UCLA_CNP", 
                          comparison_group = "ADHD", 
                          data_path = "/Users/abry4213/data/TS_feature_manuscript/",
                          plot_path = "/Users/abry4213/github/fMRI_FeaturesDisorders/plots/Manuscript_Draft/pairwise_results/",
                          node_to_from_structure = node_to_from_structure,
                            hyper_threshold = 0.6, 
                            hypo_threshold = -0.6, 
                            hyper_cmap = "Reds", 
                            hypo_cmap = "Blues")

# ABIDE ASD
plot_hyper_hypo_chord_per_group(study = "ABIDE_ASD", 
                          comparison_group = "ASD", 
                          data_path = "/Users/abry4213/data/TS_feature_manuscript/",
                          plot_path = "/Users/abry4213/github/fMRI_FeaturesDisorders/plots/Manuscript_Draft/pairwise_results/",
                          node_to_from_structure = node_to_from_structure,
                            hyper_threshold = 0.4, 
                            hypo_threshold = -0.4, 
                            hyper_cmap = "Reds", 
                            hypo_cmap = "Blues")


###############################################################################
# Glass brain visualizations
from nilearn.maskers import NiftiMasker
from nilearn import datasets, plotting
from nilearn.connectome import ConnectivityMeasure
from nilearn.maskers import MultiNiftiLabelsMasker

# Extract coordinates form Yeo atlas parcellations
# grab center coordinates for atlas labels
coordinates_dk = plotting.find_parcellation_cut_coords(labels_img="/Users/abry4213/data/neuroimaging_atlases/mni152_space/atlas-desikankilliany.nii.gz",
                                                       return_label_names=False,
                                                       background_label=0)[0:82]
np.savetxt("/Users/abry4213/data/TS_feature_manuscript/aparc_aseg_MNI_coords.csv", coordinates_dk, delimiter=",")

# Harvard-Oxford atlas
dataset_ho = datasets.fetch_atlas_harvard_oxford("cort-maxprob-thr25-1mm")
coordinates_ho = plotting.find_parcellation_cut_coords(labels_img=dataset_ho["maps"])
np.savetxt("/Users/abry4213/data/TS_feature_manuscript/HO_MNI_coords.csv", coordinates_ho, delimiter=",")



### SCZ hyper-connectivity
def plot_glass_brain_by_group(study, comparison_group):
    node_colors = pd.read_csv(f"/Users/abry4213/data/TS_feature_manuscript/{study}_{comparison_group}_node_colors.csv")
    hyper_matrix = pd.read_csv(f"/Users/abry4213/data/TS_feature_manuscript/{study}_{comparison_group}_hyper_FC.csv", header=None).to_numpy()
    hyper_coords = pd.read_csv(f"/Users/abry4213/data/TS_feature_manuscript/{study}_{comparison_group}_hyper_FC_coords.csv", header=None).to_numpy()
    hypo_matrix = pd.read_csv(f"/Users/abry4213/data/TS_feature_manuscript/{study}_{comparison_group}_hypo_FC.csv", header=None).to_numpy()
    hypo_coords = pd.read_csv(f"/Users/abry4213/data/TS_feature_manuscript/{study}_{comparison_group}_hypo_FC_coords.csv", header=None).to_numpy()
    
    
    # Plot SCZ DK connectivity
    plotting.plot_connectome(
        hyper_matrix,
        hyper_coords,
        output_file = f"/Users/abry4213/github/fMRI_FeaturesDisorders/plots/Manuscript_Draft/pairwise_results/{study}_{comparison_group}_hyper_brain.svg",
        display_mode="xz",
        node_size=25,
        node_color = node_colors.Hyper_Nodes.to_list()
    )
    plotting.plot_connectome(
        hypo_matrix,
        hypo_coords,
        output_file = f"/Users/abry4213/github/fMRI_FeaturesDisorders/plots/Manuscript_Draft/pairwise_results/{study}_{comparison_group}_hypo_brain.svg",
        display_mode="xz",
        node_size=25,
        node_color = node_colors.Hypo_Nodes.to_list()
    )
    
plot_glass_brain_by_group("UCLA_CNP", "SCZ")
plot_glass_brain_by_group("UCLA_CNP", "BPD")
plot_glass_brain_by_group("UCLA_CNP", "ADHD")
plot_glass_brain_by_group("ABIDE_ASD", "ASD")
