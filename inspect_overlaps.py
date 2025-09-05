#!/usr/bin/env python
# coding: utf-8

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from skimage.io import imread
import napari
import xarray as xr
from starfish import Experiment
from starfish.core.image._registration.transforms_list import TransformsList
from starfish.image import ApplyTransform
from starfish.types import Axes
import argparse
from matplotlib.colors import ListedColormap
from magicgui import magicgui
import json
from starfish.types import Axes, Levels
from starfish.image import ApplyTransform, LearnTransform, Filter
from starfish import ImageStack
import json

def point_color_from_target(target):
    if "missing" in str(target):
        return "gray"
    elif "invalid" in str(target):
        return "red"
    else:
        return "white"

def common_roi(transforms_list):
    # crop the images to the common region
    translations = np.array([x[2].translation for x in transforms_list.transforms])
    lower = -np.floor(translations.min(axis=0)).astype(int)
    upper = -np.ceil(translations.max(axis=0)).astype(int)
    # the `or None` is to handle the case where no cropping is necessary at the upper bound (see: https://stackoverflow.com/a/11337953/4969760)
    return slice(lower[0],upper[0] or None), slice(lower[1],upper[1] or None)

parser = argparse.ArgumentParser(description='Inspect results of starfish pipeline')
parser.add_argument('--hpi', type=int, default=0, help='hours post infection')
parser.add_argument('--fov1', type=int, default=0, help='index of the field of view 1')
parser.add_argument('--fov2', type=int, default=1, help='index of the field of view 2')
parser.add_argument('--rep', type=int, default=0, help='replication')
args = parser.parse_args()

rep = args.rep
hpi = args.hpi
fov1 = args.fov1
fov2 = args.fov2

data_folder = f"/tank/s391697/in-situ-seq-influenza/data/starfish/3nt/PR8/rep{rep}/{hpi}hpi"
prefix1 = f"pipeline_output/final_pipe/3nt/PR8/rep{rep}/{hpi}hpi/fov_{fov1}"
prefix2 = f"pipeline_output/final_pipe/3nt/PR8/rep{rep}/{hpi}hpi/fov_{fov2}"
ov_coords = pd.read_csv("data/overlap_coords.csv").query("rep ==@rep & hpi==@hpi & fov1==@fov1 & fov2==@fov2").iloc[0]
spots = pd.read_csv("data/all_spots.csv", keep_default_na=False, low_memory=False).query("rep==@rep & hpi==@hpi")
spots1 = spots.query("fov==@fov1 & x>=@ov_coords.x_min_px_fov1 & x<=@ov_coords.x_max_px_fov1 & y>=@ov_coords.y_min_px_fov1 & y<=@ov_coords.y_max_px_fov1").copy()
spots2 = spots.query("fov==@fov2 & x>=@ov_coords.x_min_px_fov2 & x<=@ov_coords.x_max_px_fov2 & y>=@ov_coords.y_min_px_fov2 & y<=@ov_coords.y_max_px_fov2").copy()
if len(spots1) > 0:
    spots1["z"] = 0
if len(spots2) > 0:
    spots2["z"] = 0

viewer = napari.Viewer()

## load images
primary1 = xr.open_dataset(f"{prefix1}_primary.nc").__xarray_dataarray_variable__
primary2 = xr.open_dataset(f"{prefix2}_primary.nc").__xarray_dataarray_variable__

y0_fov1 = ov_coords.y_min_px_fov1
x0_fov1 = ov_coords.x_min_px_fov1
y0_fov2 = ov_coords.y_min_px_fov2
x0_fov2 = ov_coords.x_min_px_fov2
dy = y0_fov1 - y0_fov2
dx = x0_fov1 - x0_fov2
translate2 = (0, dy, dx)

for i,b,c in zip(range(3,-1,-1), "AGTC"[::-1], ["magenta", "green", "yellow", "red"][::-1]):
    viewer.add_image(primary1[:,i], name=b + f" fov{fov1}", colormap=c, blending="additive", contrast_limits=[0, 0.4], visible=True)

for i,b,c in zip(range(3,-1,-1), "AGTC"[::-1], ["magenta", "green", "yellow", "red"][::-1]):
    viewer.add_image(primary2[:,i], name=b + f" fov{fov2}", colormap=c, blending="additive", contrast_limits=[0, 0.4], visible=True, translate=translate2)

## overlap rectangle
rect = [[(ov_coords.y_min_px_fov1, ov_coords.x_min_px_fov1), (ov_coords.y_max_px_fov1, ov_coords.x_max_px_fov1)]]
viewer.add_shapes(
    rect,
    shape_type="rectangle",
    edge_color="red",
    face_color="transparent",
    opacity=0.5,
    name=f"Overlap FOV {fov1}/{fov2}",
    visible=False
)

## Registration overlaps
crop_img1 = primary1.isel(  ## crop overlaps
    y=slice(int(round(ov_coords.y_min_px_fov1)), int(round(ov_coords.y_max_px_fov1))),
    x=slice(int(round(ov_coords.x_min_px_fov1)), int(round(ov_coords.x_max_px_fov1)))
)
crop_img2 = primary2.isel(
    y=slice(int(round(ov_coords.y_min_px_fov2)), int(round(ov_coords.y_max_px_fov2))),
    x=slice(int(round(ov_coords.x_min_px_fov2)), int(round(ov_coords.x_max_px_fov2)))
)

max1 = ImageStack.from_numpy(crop_img1)  ## MIP of overlps
max2= ImageStack.from_numpy(crop_img2)
max1 = max1.reduce({Axes.CH, Axes.ROUND}, func="max").xarray.squeeze()
max2 = max2.reduce({Axes.CH, Axes.ROUND}, func="max").xarray.squeeze()

y_common, x_common = min(max1.shape[0], max2.shape[0]), min(max1.shape[1], max2.shape[1])
crop_max1 = max1.isel(x=slice(0, x_common), y=slice(0, y_common))
crop_max2 = max2.isel(x=slice(0, x_common), y=slice(0, y_common))

stack_np = np.stack([crop_max1.values, crop_max2.values], axis=0)
stack_np = stack_np[:, np.newaxis, np.newaxis, :, :]
stack_both = ImageStack.from_numpy(stack_np)

learn_translation = LearnTransform.Translation(reference_stack=stack_both.sel({Axes.ROUND: 0}), axes=Axes.ROUND, upsampling=1000) ## Registration
warp = ApplyTransform.Warp()
transform_list = learn_translation.run(stack_both)
registered_stack = warp.run(stack_both, transforms_list=transform_list)
transform_list.to_json(f"test_reg.json")

with open(f"test_reg.json", "r", encoding="utf-8") as file_reg:
    dict_reg = json.load(file_reg)
tx = dict_reg["transforms_list"][1][2][0][2]
ty = dict_reg["transforms_list"][1][2][1][2]
tx *= -1
ty *= -1

## add registered MIP overlaps
viewer.add_image(crop_max1, blending="additive", colormap="red", name="max1", visible=False, translate=(ov_coords.y_min_px_fov1, ov_coords.x_min_px_fov1))
viewer.add_image(crop_max2, blending="additive", colormap="red", name="max2", visible=False, translate=(ov_coords.y_min_px_fov1, ov_coords.x_min_px_fov1))
viewer.add_image(crop_max2, blending="additive", colormap="red", name="max 2 trans", visible=False, translate=(ov_coords.y_min_px_fov1 + ty, ov_coords.x_min_px_fov1 + tx))

## add MIP outlines
crop_max2_trans = crop_max2.assign_coords(
    x = crop_max2.coords["x"] + tx,
    y = crop_max2.coords["y"] + ty,
)

for img, name in zip([crop_max1, crop_max2_trans], ["rect max1", "rect max2 trans"]):
    rect_max = [[(min(img.coords["y"]), min(img.coords["x"])), (max(img.coords["y"]), max(img.coords["x"]))]]
    viewer.add_shapes(
        rect_max,
        shape_type="rectangle",
        edge_color="green",
        face_color="transparent",
        opacity=0.5,
        name=name,
        translate=(ov_coords.y_min_px_fov1, ov_coords.x_min_px_fov1),
        edge_width=5,
        visible=False
    )

## Common overlap
h, w = crop_max1.shape

y0 = max(0, 0)         
y1 = min(h, h)
y0 = max(y0, ty)
y1 = min(y1, h+ty)

x0 = max(0, 0)
x1 = min(w, w)
x0 = max(x0, tx)
x1 = min(x1, w+tx)

reg1 = registered_stack.xarray.sel({Axes.ROUND: 0}).squeeze().values
reg2 = registered_stack.xarray.sel({Axes.ROUND: 1}).squeeze().values
reg1 = reg1.astype(float)
reg2 = reg2.astype(float)

common1 = reg1[int(y0):int(y1), int(x0):int(x1)]
common2 = reg2[int(y0):int(y1), int(x0):int(x1)]

g_y0 = ov_coords.y_min_px_fov1 + y0
g_x0 = ov_coords.x_min_px_fov1 + x0
g_y1 = ov_coords.y_min_px_fov1 + y1
g_x1 = ov_coords.x_min_px_fov1 + x1

# Rechteck
rect = [[(g_y0, g_x0), (g_y1, g_x1)]]
viewer.add_shapes(
    rect,
    shape_type="rectangle",
    edge_color="blue",
    edge_width=5,
    face_color="transparent",
    name="common overlap",
    visible=True
)

# beide Overlap‐MIPs einblenden
viewer.add_image(
    common1,
    name="common1",
    colormap="red",
    blending="additive",
    translate=(g_y0, g_x0),#(ov_coords.y_min_px_fov1 + y0, ov_coords.x_min_px_fov1 + x0),
    visible=False
)
viewer.add_image(
    common2,
    name="common2 registered",
    colormap="red",
    blending="additive",
    translate=(g_y0, g_x0),#(ov_coords.y_min_px_fov1 + y0, ov_coords.x_min_px_fov1 + x0),
    visible=False
)

# load spots
spots1_ov = spots1.query("x>=@g_x0 & x<=@g_x1 & y>=@g_y0 & y<=@g_y1").copy()
spots2_ov = spots2.copy()
spots2_ov["x"] = spots2_ov["x"] + (dx + tx)
spots2_ov["y"] = spots2_ov["y"] + (dy + ty)
spots2_ov = spots2_ov.query("x>=@g_x0 & x<=@g_x1 & y>=@g_y0 & y<=@g_y1")

colors = ["red", "blue"]
for i, (spots_df, fov_id) in enumerate(zip([spots1_ov, spots2_ov], [fov1, fov2])):
    spots_np = spots_df[['z','y','x','radius']].astype(float).to_numpy()
    features = {'target': spots_df['target'].replace("missing","").values}

    viewer.add_points(
        spots_np[:, :3],  
        size=spots_np[:, 3],
        name=f"Spots FOV {fov_id}",
        symbol="ring",
        features=features,
        text="target",
        face_color=colors[i],
        visible=False,
        border_width=0.5
    )

## matched spots
overlap = f"{fov1}/{fov2}"
spots_match = pd.read_csv("data/overlap_matches_with_distances_sbd.csv", keep_default_na=False).query("rep==@rep & hpi==@hpi & overlap==@overlap")
spots_match["z"] = 0
spots_match["radius"] = 3

matched1 = spots_match.query("match==True & fov1==@fov1")[['z','y_fov1','x_fov1','radius']].to_numpy()
matched2 = spots_match.query("match==True & fov2==@fov2")[['z','y_fov2','x_fov2','radius']].to_numpy()
unmatched1 = spots_match.query("match==False & x_fov2==-1")[['z','y_fov1','x_fov1','radius']].to_numpy()
unmatched2 = spots_match.query("match==False & x_fov1==-1")[['z','y_fov2','x_fov2','radius']].to_numpy()

for df, name in zip([matched1, matched2], [f"matched spots fov{fov1}", f"matched spots fov{fov2}"]):
    viewer.add_points(
        df[:, :3],  
        size=df[:, 3],
        name=name,
        symbol="ring",
        face_color="green",
        edge_color="green",
        visible=True,
        border_width=0.5
    )

for df, name in zip([unmatched1, unmatched2], [f"unmatched spots fov{fov1}", f"unmatched spots fov{fov2}"]):
    viewer.add_points(
        df[:, :3],  
        size=df[:, 3],
        name=name,
        symbol="ring",
        face_color="red",
        edge_color="red",
        visible=True,
        border_width=0.5
    )

napari.run()