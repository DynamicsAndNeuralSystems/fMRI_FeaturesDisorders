import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import nxviz as nv

# Create a brain region interaction df
df_brain = pd.read_csv("/Users/abry4213/data/TS_feature_manuscript/UCLA_CNP_SCZ_hyper_FC_pairs.csv")
node_to_from_structure =  pd.read_csv("/Users/abry4213/data/TS_feature_manuscript/node_to_from_structure.csv")
interactions_brain = np.array(df_brain)
G_brain = nx.DiGraph(name = "Schizophrenia Hyperconnectivity")
for i in range(len(interactions_brain)):
    interaction = interactions_brain[i]
    a = interaction[0] # protein a node
    b = interaction[1] # protein b node
    w = interaction[2] # score as weighted edge
    
    # To only keep high scoring edges, use the following lines
    G_brain.add_weighted_edges_from([(a,b,w)])
    
for n, d in G_brain.nodes(data=True):
    brain_region = n
    lobe = node_to_from_structure.query("to == @brain_region")["from"].to_string(index=False)
    G_brain.nodes[n]["lobe"] = lobe
    
for u, v, d in G_brain.edges(data=True):
    print(u)
    print(v)
    print(d)
    # d["weight"] = randint(0, 10)
    
    
nv.circos(G_brain, group_by="lobe", node_color_by="lobe", 
          edge_color_by="weight", edge_width=1)


########################

from nxviz import lines, utils, layouts, nodes, plots, edges
from random import choice

G = nx.erdos_renyi_graph(n=71, p=0.1)
for n, d in G.nodes(data=True):
    G.nodes[n]["group"] = choice(["a", "b", "c"])
    G.nodes[n]["value"] = np.random.exponential()
    

for u, v, d in G.edges(data=True):
    G.edges[u, v]["edge_value"] = np.random.exponential()
    
# Customize node styling
ax = plt.gca()

nt = utils.node_table(G)
pos = layouts.circos(nt, group_by="group", sort_by="value")

from nxviz import encodings as aes


def group_colormap(data: pd.Series):
    cmap = {"a": "black", "b": "blue", "c": "red"}
    return data.apply(lambda x: cmap.get(x))


def value_colormap(data: pd.Series):
    """Value colormap."""
    norm = plt.cm.Normalize(vmin=data.min(), vmax=data.max())
    cmap = plt.cm.get_cmap("viridis")
    return data.apply(lambda x: cmap(norm(x)))


def node_size(data: pd.Series):
    return data.apply(np.sqrt)

node_color = group_colormap(nt["group"])
alpha = nodes.transparency(nt, alpha_by=None)
size = nodes.node_size(nt, "value")
patches = nodes.node_glyphs(
    nt, pos, node_color=node_color, alpha=alpha, size=size
)
for patch in patches:
    ax.add_patch(patch)

# Customize edge styling
et = utils.edge_table(G)
edge_color = edges.edge_colors(et, nt=None, color_by=None, node_color_by=None)
lw = edges.line_width(et, lw_by=None)
alpha = edges.transparency(et, alpha_by="edge_value")
patches = lines.circos(
    et, pos, edge_color=edge_color, alpha=alpha, lw=lw, aes_kw={"fc": "none"}
)
for patch in patches:
    ax.add_patch(patch)

plots.rescale(G)
plots.aspect_equal()
plots.despine()