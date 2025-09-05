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


parser = argparse.ArgumentParser(description='Inspect results of starfish pipeline')
parser.add_argument('--hpi', type=int, default=5, help='hours post infection')
parser.add_argument('--fov_index', type=int, default=0, help='index of the field of view to inspect')
parser.add_argument('--strain', type=str, default="PR8", help='virus strain')
parser.add_argument('--barcode', type=str, default="3nt", help='barcode used in the experiment')
parser.add_argument('--rep', type=int, default=0, help='replication')
parser.add_argument('--nucleus-segmentation-suffix', type=str, default="_nuclei.png", help='file name suffix of the nucleus segmentation mask')
parser.add_argument('--cell-segmentation-suffix', type=str, default="_cells.png", help='file name suffix of the cell segmentation mask')
parser.add_argument('--exp_name', type=str, help='name of experiment', default="final_pipe")
parser.add_argument("--decoding", type=str, help="decoding strategie", default="_sbd")
args = parser.parse_args()

hpi = args.hpi
fov_index = args.fov_index
strain = args.strain
barcode = args.barcode
rep = args.rep
decoding_strategie = f"{args.decoding}"
nucleus_segmentation_suffix = args.nucleus_segmentation_suffix
cell_segmentation_suffix = args.cell_segmentation_suffix
exp_name = args.exp_name
data_folder = f"/tank/s391697/in-situ-seq-influenza/data/starfish/{barcode}/{strain}/rep{rep}/{hpi}hpi"
analysis_folder = f"/tank/s391697/in-situ-seq-influenza/analysis/starfish/{barcode}/{strain}/rep{rep}/{hpi}hpi" if exp_name == "None" else f"/home/witte/master/pipeline_output/{exp_name}/{barcode}/{strain}/rep{rep}/{hpi}hpi"
prefix = f"{analysis_folder}/fov_{fov_index}"

# ## Load raw data
exp = Experiment.from_json(f"{data_folder}/experiment.json")

fov = exp.fovs()[fov_index]
imgs = fov.get_image("primary")
nuclei = fov.get_image("nuclei")
bf = fov.get_image("brightfield")


# ### Apply registration
transforms_list = TransformsList()
transforms_list = transforms_list.from_json(f'{prefix}_registration.json')

warp = ApplyTransform.Warp()
registered_imgs = warp.run(imgs, transforms_list=transforms_list)
transforms_list = transforms_list.from_json(f'{prefix}_registration.json')
registered_nuclei = warp.run(nuclei, transforms_list=transforms_list)
transforms_list = transforms_list.from_json(f'{prefix}_registration.json')
registered_bf = warp.run(bf, transforms_list=transforms_list)


# #### Crop to common region
def common_roi(transforms_list):
    # crop the images to the common region
    translations = np.array([x[2].translation for x in transforms_list.transforms])
    lower = -np.floor(translations.min(axis=0)).astype(int)
    upper = -np.ceil(translations.max(axis=0)).astype(int)
    # the `or None` is to handle the case where no cropping is necessary at the upper bound (see: https://stackoverflow.com/a/11337953/4969760)
    return slice(lower[0],upper[0] or None), slice(lower[1],upper[1] or None)


xroi, yroi = common_roi(transforms_list)
registered_imgs = registered_imgs.sel({Axes.X: xroi, Axes.Y: yroi})
registered_nuclei = registered_nuclei.sel({Axes.X: xroi, Axes.Y: yroi})
registered_bf = registered_bf.sel({Axes.X: xroi, Axes.Y: yroi})


# ## Load segmentation results
nucleus_mask = imread(f"{prefix}{nucleus_segmentation_suffix}")
cell_mask = imread(f"{prefix}{cell_segmentation_suffix}")

## If the segmentation was done on the full image, crop it to the common region
if registered_nuclei.shape[Axes.Y] != nucleus_mask.shape[0] or registered_nuclei.shape[Axes.X] != nucleus_mask.shape[1]:
    nucleus_mask = nucleus_mask[yroi, xroi]
    cell_mask = cell_mask[yroi, xroi]

# ## Load processed primary images
primary = xr.open_dataset(f"{prefix}_primary.nc").__xarray_dataarray_variable__

# ## Load detected spots
spots = pd.read_csv(f"{prefix}_spots_sbd.csv", index_col=0, keep_default_na=False)
spots_np = spots[['z', 'y','x','radius']].to_numpy()

spots_np_pixel = None
if f"fov_{fov_index}_spots_pbd.csv" in os.listdir(analysis_folder):
    spots_pixel_based = pd.read_csv(f"pipeline_output/final_pipe/3nt/PR8/rep{rep}/{hpi}hpi/fov_{fov_index}_spots_pbd.csv")
    spots_pixel_based["diameter"] = [x * 2 for x in spots_pixel_based.radius]
    spots_np_pixel = spots_pixel_based[['y','x','diameter']].to_numpy()

# ## View results in napari
def point_color_from_target(target):
    if "missing" in str(target):
        return "gray"
    elif "invalid" in str(target):
        return "red"
    else:
        return "white"

def point_color_from_iter(iter):
    colors_pbd = ["purple", "red", "orange", "yellow", "green", "blue"]
    if "," in iter:
        return "white"
    else:
        return colors_pbd[int(iter)]

def get_border_color(segment, default_color):
    """
    Set border color to red if the segment contains 'invalid' or 'missing'.
    Otherwise, return the default color.
    """
    if isinstance(segment, float) and np.isnan(segment):
        return default_color  # Behalte die ursprüngliche Farbe, wenn der Wert NaN ist
    segment_str = str(segment).strip().lower()  # Konvertiere zu String, entferne Leerzeichen und mache klein
    if "invalid" in segment_str or "missing" in segment_str:
        return [1, 0, 0, 1]  # Rot
    return default_color  # Behalte die ursprüngliche Farbe

contrast_limits = [0, 0.05]
if rep != 0:
    contrast_limits = [0, 0.002]

viewer = napari.Viewer()
for i,b,c in zip(range(3,-1,-1), "AGTC"[::-1], ["magenta", "green", "yellow", "red"][::-1]):
    viewer.add_image(registered_imgs.xarray[:,i], name=b+" (raw)", colormap=c, blending="additive", visible=False, contrast_limits=contrast_limits)
viewer.add_image(registered_bf.xarray[:,0], name="brightfield", opacity=0.2, visible=False)
viewer.add_image(registered_nuclei.xarray[:,0], name="dapi", opacity=0.2, visible=False)
for i,b,c in zip(range(3,-1,-1), "AGTC"[::-1], ["magenta", "green", "yellow", "red"][::-1]):
    viewer.add_image(primary[:,i], name=b, colormap=c, blending="additive", contrast_limits=[0, 0.4], visible=True)


## spot based 
viewer.add_points(
    spots_np[:,:3],
    size=spots_np[:,3],
    name="spots spot-based",
    symbol="ring",
    features={"target": spots.target.replace("missing", "").values},
    text="target",
    face_color=[point_color_from_target(x) for x in spots.target.values],
    visible=False
)
## pixel based
if spots_np_pixel is not None:
    viewer.add_points(
        spots_np_pixel[:,:2],
        size=spots_np_pixel[:,2],
        name="spots pixel-based",
        symbol="ring",
        features={"target": spots_pixel_based.target.replace("missing", "").values},

        text=[
            f"{t} (r={r:.2f})"
            for t, r in zip(spots_pixel_based.target.values, spots_pixel_based.radius.values)
        ],
        face_color="white",
        visible=False
    )

rect_coords = [[800, 800], [800, 950], [650, 950], [650, 800]]
viewer.add_shapes(
    [rect_coords],
    shape_type="polygon",
    edge_color="red",
    face_color="transparent",
    opacity=0.5,
    name="Annotation Area"
)

if f"rep{rep}_hpi{hpi}_fov{fov_index}.csv" in os.listdir("manual_annotation/"):
## Manuell annotierte spots
    spots_annot = pd.read_csv(f"manual_annotation/rep{rep}_hpi{hpi}_fov{fov_index}.csv", keep_default_na=False)
    if len(spots_annot) > 0:
        colors = []
        for spot in spots_annot.text:
            if "0" in spot or "1" in spot:
                colors.append("red")
            else:
                colors.append("blue")
        spots_annot_np = spots_annot[["x", "y"]].to_numpy()
        viewer.add_points(
            spots_annot_np,
            size=8,
            name="manual annotation",
            border_width=0.05,
            border_color=colors,  # Blau 
            face_color=[0, 0, 0, 0], 
        )

    ## Vergleich manuelle Annotation und Pipeline
    file_name = f"rep{rep}_hpi{hpi}_fov{fov_index}.csv"
    spots_matching_df = pd.read_csv(f"manual_annotation_decoded/matched_spots_{exp_name}{decoding_strategie}.csv")
    spots_matching_df = spots_matching_df[spots_matching_df["file"] == file_name]

    matched_spots = spots_matching_df[~ spots_matching_df["distance"].isna()]
    unmatched_spots_pipe = spots_matching_df[spots_matching_df["x_man"].isna()]
    unmatched_spots_man = spots_matching_df[spots_matching_df["x_pipe"].isna()]

    unique_colors = np.random.rand(len(matched_spots),3)
    cmap = ListedColormap(unique_colors)
    values = np.arange(len(matched_spots))
    colors = cmap(values)

    if len(matched_spots) > 0:
        viewer.add_points( #manuelle spots
            matched_spots[["y_man", "x_man"]].to_numpy(),
            size = 5,
            name="manual annotated spots matched",
            border_width=0.05,
            border_color=[get_border_color(seg, c) for seg, c in zip(matched_spots["segments_man"], colors)],
            face_color=[0,0,0,0]
        )

        viewer.add_points( #pipeline spots
            matched_spots[["y_pipe", "x_pipe"]].to_numpy(),
            size = 5,
            name="Pipeline annotated spots matched",
            border_width=0.05,
            border_color=[get_border_color(seg, c) for seg, c in zip(matched_spots["segments_pipe"], colors)],
            face_color=[0,0,0,0]
        )

    if len(unmatched_spots_pipe) > 0:
        viewer.add_points( #pipeline unmatched spots
            unmatched_spots_pipe[["y_pipe", "x_pipe"]].to_numpy(),
            size = 5,
            name="Pipeline unmatched spots",
            border_width=0.05,
            border_color=[get_border_color(seg, [0.6, 0.6, 0.6]) for seg in unmatched_spots_pipe["segments_pipe"]],
            face_color=[0,0,0,0]
        )

    if len(unmatched_spots_man) > 0:
        viewer.add_points( #manuell unmatched spots
            unmatched_spots_man[["y_man", "x_man"]].to_numpy(),
            size = 5,
            name="Manual unmatched spots",
            border_width=0.05,
            border_color=[get_border_color(seg, [1.0, 1.0, 1.0]) for seg in unmatched_spots_pipe["segments_man"]],
            face_color=[0,0,0,0]
        )

## density mask
density_map = np.load(f"pipeline_output/mask_percentile_10/{barcode}/{strain}/rep{rep}/{hpi}hpi/fov_{fov_index}_density_map.npy")

viewer.add_image(density_map, name="Density map", colormap='gray', visible=False)

# prepare empty mask layer
mask_layer = viewer.add_labels(np.zeros_like(density_map, dtype=int), name="Mask", opacity=0.4, colormap={1: 'white'})

@magicgui(
    percentile={"widget_type": "FloatSlider", "min": 0.0, "max": 100.0, "step": 1.0}
)
def update_mask(percentile: float = 50.0):
    threshold_value = np.percentile(density_map, percentile)
    mask = (density_map >= threshold_value).astype(int)
    mask_layer.data = mask

viewer.window.add_dock_widget(update_mask)


path = f"pipeline_output/masks_threshold_model/{barcode}/{strain}/rep{rep}/{hpi}hpi/fov_{fov_index}_density_mask"

mask = np.load(f"{path}.npy")
viewer.add_image(mask, name=f"density_mask", visible=False, opacity=0.4)
mask_filtered = np.load(f"{path}_filtered.npy")
viewer.add_image(mask_filtered, name="density_mask_filtered", visible=False, opacity=0.4, colormap="red")

# calculate area for zoom
min_coords = np.min(rect_coords, axis=0)
max_coords = np.max(rect_coords, axis=0)
center_coords = (min_coords + max_coords) / 2 
zoom_level = 2.0

# focus camera on rectangle
viewer.camera.center = center_coords  
viewer.camera.zoom = zoom_level     

viewer.add_labels(nucleus_mask, name="nuclei", opacity=0.4, visible=False)
viewer.add_labels(cell_mask, name="cells", opacity=0.2, visible=False)

mip = primary.max(dim=["r", "c","z"])
viewer.add_image(mip, visible=False)

napari.run()