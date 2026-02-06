import argparse
import napari
import argparse
import xarray as xr
from starfish.types import Axes
from starfish.core.imagestack.parser.numpy import NumpyData
from starfish.core.imagestack.imagestack import ImageStack
import numpy as np
from starfish.image import ApplyTransform, LearnTransform, Filter
import json
import glob
import pandas as pd
import networkx as nx
from magicgui import magicgui
from napari.layers import Image


## Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--strain", default="PR8", type=str)
parser.add_argument("--rep", default=0, type=int)
parser.add_argument("--hpi", default=5, type=int)
parser.add_argument("--python_version", type=str)
args = parser.parse_args()

## Variables
strain = args.strain
rep = args.rep 
hpi = args.hpi
python_version = args.python_version
base_dir = "/data/influenza-genome-packaging/results/preprocessed"

## Functions
from functions import get_range, overlap_1d, overlap_2d, registration, crop_overlap

## Load images and MIP
fov_paths = sorted(glob.glob(f"{base_dir}/{strain}/rep{rep}/{hpi}hpi/*"))
fovs = []
mps = []

for path in fov_paths:
    fov = xr.open_dataset(path).__xarray_dataarray_variable__
    mp = fov.max(dim=["r", "c"]).squeeze()
    fovs.append(fov)
    mps.append(mp)

## Calculate pixel size
px_size_x = fov[0].xc.values[1] - fov[0].xc.values[0]
px_size_y = fov[0].yc.values[1] - fov[0].yc.values[0]

for fov in fovs:
    if fov.xc.values[1] - fov.xc.values[0] != px_size_x:
        print("WARNING: x pixel size differs between fovs")
    
    if fov.yc.values[1] - fov.yc.values[0] != px_size_y:
        print("WARNING: y pixel size differs between fovs")

## Add FOVs
viewer = napari.Viewer()

colors = ["red", "blue", "green", "yellow", "pink", "magenta", "brown", "gray", "cyan", "lime", "purple", "orange", "red", "blue", "green", "yellow", "pink", "magenta", "brown", "gray", "cyan", "lime", "purple", "orange"]
scale = [px_size_y, px_size_x]
translates = []

for mp, mp_idx, c in zip(mps, range(0, len(mps)), colors[:len(mps)]):
    x0 = mp.xc.values[0]
    y0 = mp.yc.values[0]
    translate = [y0, x0]
    translates.append(translate)

    viewer.add_image(mp, name=f"mp{mp_idx}", colormap=c, blending="additive", contrast_limits=[0, 0.4], visible=False, translate=translate, scale=scale)

## Load MST and registration data
regs = pd.read_csv(f"edges_and_registration_py{python_version}.csv").query("strain==@strain & rep==@rep & hpi==@hpi")
edges = [(source, target, weight) for source, target, weight in zip(list(regs["source"]), list(regs["target"]), list(regs["weight"]))]

# Identify connected MSTs and root for each MST
G = nx.Graph()
G.add_weighted_edges_from(edges)

for mst in nx.connected_components(G):
    sub = G.subgraph(mst)

    node_weight_sum = {
        n: sum(d.get("weight", 0) for _, _, d in sub.edges(n, data=True))
        for n in sub.nodes
    }

    root = list(node_weight_sum.keys())[list(node_weight_sum.values()).index(max(node_weight_sum.values()))]
    viewer.add_image(mps[root], name=f" mp{root} reg", colormap=colors[root], blending="additive", contrast_limits=[0, 0.4], visible=True, translate=translates[root], scale=scale)
    tree = nx.bfs_tree(sub, root)

    branches = {
        child: list(nx.descendants(tree, child)) + [child]
        for child in tree.successors(root)
    }

    # print(regs)
    # print("root: ", root)

    for branch in branches.items():
        # print("branch: ", branch)
        nodes = branch[1][::-1]
        prev_node = root
        tx = 0
        ty = 0

        for node in nodes:
            # print("node: ", node)
            # print("prev_node: ", prev_node)

            source_target_df1 = regs.query("source==@prev_node & target==@node")
            source_target_df2 = regs.query("source==@node & target==@prev_node")

            if len(source_target_df1) == 1:
                tx += source_target_df1["tx"].values[0]
                ty += source_target_df1["ty"].values[0]

            elif len(source_target_df2) == 1:
                tx += (source_target_df2["tx"].values[0] * -1)
                ty += (source_target_df2["ty"].values[0] * -1)

            # print("df1: ", source_target_df1)
            # print("df2: ", source_target_df2)
            # print("Registration: ", tx, ty)

            trans = translates[node]
            trans_y = trans[0] + ty
            trans_x = trans[1] + tx

            # print("Translate: ", trans)

            viewer.add_image(mps[node], name=f" mp{node} reg", colormap=colors[node], blending="additive", contrast_limits=[0, 0.4], visible=True, translate=[trans_y, trans_x], scale=scale)
            prev_node = node

@magicgui(call_button="Toggle Registration")
def toggle_images():
    for layer in viewer.layers:
        if isinstance(layer, Image) and "reg" in layer.name:
            layer.visible = not layer.visible

        if isinstance(layer, Image) and "reg" not in layer.name:
            layer.visible = not layer.visible

viewer.window.add_dock_widget(toggle_images)
napari.run()