# import argparse
# import napari
# import argparse
# import xarray as xr
# from starfish.types import Axes
# from starfish.core.imagestack.parser.numpy import NumpyData
# from starfish.core.imagestack.imagestack import ImageStack
# import numpy as np
# from starfish.image import ApplyTransform, LearnTransform, Filter
# import json


# ## Arguments
# parser = argparse.ArgumentParser()
# parser.add_argument("--strain", default="PR8", type=str)
# parser.add_argument("--rep", default=0, type=int)
# parser.add_argument("--hpi", default=5, type=int)
# parser.add_argument("--fov1", default=0, type=int)
# parser.add_argument("--fov2", default=1, type=int)
# args = parser.parse_args()

# ## Variables
# strain = args.strain
# rep = args.rep 
# hpi = args.hpi
# fov1_idx = args.fov1
# fov2_idx = args.fov2
# base_dir = "/data/influenza-genome-packaging/results/preprocessed"
# path_fov1 = f"{base_dir}/{strain}/rep{rep}/{hpi}hpi/primary-fov_00{fov1_idx}.nc" if fov1_idx < 10 else f"{base_dir}/{strain}/rep{rep}/{hpi}hpi/primary-fov_0{fov1_idx}.nc"
# path_fov2 = f"{base_dir}/{strain}/rep{rep}/{hpi}hpi/primary-fov_00{fov2_idx}.nc" if fov2_idx < 10 else f"{base_dir}/{strain}/rep{rep}/{hpi}hpi/primary-fov_0{fov2_idx}.nc"

# ## Functions
# from functions import get_range, overlap_1d, overlap_2d, registration, crop_overlap

# ## Load images and MIP
# fov1 = xr.open_dataset(path_fov1).__xarray_dataarray_variable__
# fov2 = xr.open_dataset(path_fov2).__xarray_dataarray_variable__
# mp1 = fov1.max(dim=["r", "c"]).squeeze()
# mp2 = fov2.max(dim=["r", "c"]).squeeze()

# ## Calculate pixel size
# px_size_x1 = fov1.xc.values[1] - fov1.xc.values[0]
# px_size_y1 = fov1.yc.values[1] - fov1.yc.values[0]
# px_size_x2 = fov2.xc.values[1] - fov2.xc.values[0]
# px_size_y2 = fov2.yc.values[1] - fov2.yc.values[0]

# if px_size_x1 != px_size_x2 or px_size_y1 != px_size_y2 or px_size_x1 != px_size_x2:
#     print("WARNING: pixel sizes differ between fov1 and fov2")

# ## Add channel and round images
# viewer = napari.Viewer()

# translates = [[0, fov1.yc[0], fov1.xc[0]], [0, fov2.yc[0], fov2.xc[0]]]

# for mp, fov_idx in [(fov1, 1), (fov2, 2)]:
#     scale = [1, px_size_y1, px_size_x1]

#     for i,b,c in zip(range(3,-1,-1), "AGTC"[::-1], ["magenta", "green", "yellow", "red"][::-1]):
#         viewer.add_image(mp[:,i], name=b + f" mp{fov_idx}", colormap=c, blending="additive", contrast_limits=[0, 0.4], visible=False, translate=translates[fov_idx-1], scale=scale)

# ## Add MIPs of FOV1 and FOV2
# for mp, mp_idx, c in [(mp1, 1, "green"), (mp2, 2, "red")]:
#     x0 = mp.xc[0]
#     y0 = mp.yc[0]

#     translate = [y0, x0]
#     scale = [px_size_y1, px_size_x1]
#     viewer.add_image(mp, name=f" mp{mp_idx}", colormap=c, blending="additive", contrast_limits=[0, 0.4], visible=False, translate=translate, scale=scale)

# ## Add Overlap Rectangle
# # Calculate FOV ranges, widths and heights
# range_fov1, range_fov2 = get_range(fov1), get_range(fov2)

# width_fov1 = max(range_fov1[0]) - min(range_fov1[0])
# height_fov1 = max(range_fov1[1]) - min(range_fov1[1])
# width_fov2 = max(range_fov2[0]) - min(range_fov2[0])
# height_fov2 = max(range_fov2[1]) - min(range_fov2[1])

# # Calculate overlap range, width and height
# overlap = overlap_2d(range_fov1, range_fov2)
# x_overlap, y_overlap = overlap

# width_overlap = max(x_overlap) - min(x_overlap)
# height_overlap = max(y_overlap) - min(y_overlap)
# width_overlap_px = width_overlap / px_size_x1
# height_overlap_px = height_overlap / px_size_y1

# # Add Rectangle
# rect_ov = [[
#     (0, 0),
#     (height_overlap_px, width_overlap_px)
# ]]

# viewer.add_shapes(
#     rect_ov,
#     shape_type="rectangle",
#     edge_color="yellow",
#     edge_width=5,
#     face_color="transparent",
#     opacity=0.5,
#     name="Overlap Area",
#     scale=[px_size_y1, px_size_x1], 
#     translate=[fov1.yc.values[0] + height_fov1 - height_overlap, fov1.xc.values[0] + width_fov1 - width_overlap]
# )

# ## Registration
# # Cut out the overlap area of both MIPs
# mps_ov = []

# for i, mp, fov_idx, c in [(1, mp1, fov1_idx, "green"), (2, mp2, fov2_idx, "red")]:
#     x_idx = np.where((mp.xc.values >= x_overlap[0]) & (mp.xc.values <= x_overlap[1]))[0]
#     y_idx = np.where((mp.yc.values >= y_overlap[0]) & (mp.yc.values <= y_overlap[1]))[0]

#     fov_ov = mp.sel(
#         x = slice(x_idx[0], x_idx[-1]),
#         y = slice(y_idx[0], y_idx[-1])
#     )

#     mps_ov.append(fov_ov)

#     viewer.add_image(fov_ov, translate=[fov_ov.yc.values[0], fov_ov.xc.values[0]], scale=scale, name=f"FOV{i} overlap", colormap=c, blending="additive")

# # Registration
# mp_ov1, mp_ov2 = mps_ov
# # stack_np = np.stack([mp_ov1.values, mp_ov2.values], axis=0)
# # stack_np = stack_np[:, np.newaxis, np.newaxis, :, :]
# # stack_both = ImageStack.from_numpy(stack_np)

# # learn_translation = LearnTransform.Translation(reference_stack=stack_both.sel({Axes.ROUND: 0}), axes=Axes.ROUND, upsampling=1000) ## Registration
# # warp = ApplyTransform.Warp()
# # transform_list = learn_translation.run(stack_both)
# # registered_stack = warp.run(stack_both, transforms_list=transform_list)
# # transform_list.to_json(f"registration.json")

# # with open(f"registration.json", "r", encoding="utf-8") as file_reg:
# #     dict_reg = json.load(file_reg)

# # tx = dict_reg["transforms_list"][1][2][0][2] * -1 * px_size_x1
# # ty = dict_reg["transforms_list"][1][2][1][2] * -1 * px_size_y1

# # viewer.add_image(fov_ov, translate=[y_overlap[0] + ty, x_overlap[0] + tx], scale=scale, name=f"FOV2 registered", colormap="yellow", blending="additive")

# tx, ty = registration(mp_ov1, mp_ov2, "registration_napari.json")
# napari.run()

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
args = parser.parse_args()

## Variables
strain = args.strain
rep = args.rep 
hpi = args.hpi
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

colors = ["red", "blue", "green", "yellow", "pink", "magenta", "brown", "gray", "cyan", "lime", "purple", "orange"]
scale = [px_size_y, px_size_x]
translates = []

for mp, mp_idx, c in zip(mps, range(0, len(mps)), colors[:len(mps)]):
    x0 = mp.xc.values[0]
    y0 = mp.yc.values[0]
    translate = [y0, x0]
    translates.append(translate)

    viewer.add_image(mp, name=f"mp{mp_idx}", colormap=c, blending="additive", contrast_limits=[0, 0.4], visible=False, translate=translate, scale=scale)

## Load MST and registration data
regs = pd.read_csv("edges_and_registration.csv").query("strain==@strain & rep==@rep & hpi==@hpi")
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

    print(regs)
    print("root: ", root)

    for branch in branches.items():
        print("branch: ", branch)
        nodes = branch[1][::-1]
        prev_node = root
        tx = 0
        ty = 0

        for node in nodes:
            print("node: ", node)
            print("prev_node: ", prev_node)
            print(regs.query("source==@prev_node & target==@node"))

            source_target_df = pd.concat([regs.query("source==@prev_node & target==@node"), regs.query("source==@node & target==@prev_node")])

            print("len_df: ", len(source_target_df))

            if len(source_target_df) < 1:
                print("WARNING: more than one entries for source-target pair")

            tx += source_target_df["tx"].values[0]
            ty += source_target_df["ty"].values[0]

            print("Registration: ", tx, ty)

            trans = translates[node]
            trans_y = trans[0] + ty
            trans_x = trans[1] + tx

            print("Translate: ", trans)

            viewer.add_image(mps[node], name=f" mp{node} reg", colormap=colors[node], blending="additive", contrast_limits=[0, 0.4], visible=True, translate=[trans_y, trans_x], scale=scale)
            prev_node = node

@magicgui(call_button="Apply Registration")
def toggle_images():
    for layer in viewer.layers:
        if isinstance(layer, Image) and "reg" in layer.name:
            layer.visible = not layer.visible

        if isinstance(layer, Image) and "reg" not in layer.name:
            layer.visible = not layer.visible

viewer.window.add_dock_widget(toggle_images)
napari.run()