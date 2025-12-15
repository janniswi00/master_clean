import argparse
import napari
import argparse
import xarray as xr
from starfish.types import Axes
from starfish.core.imagestack.parser.numpy import NumpyData
from starfish.core.imagestack.imagestack import ImageStack
import numpy as np

## Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--strain", default="PR8", type=str)
parser.add_argument("--rep", default=0, type=int)
parser.add_argument("--hpi", default=5, type=int)
parser.add_argument("--fov1", default=0, type=int)
parser.add_argument("--fov2", default=1, type=int)
args = parser.parse_args()

## Variables
strain = args.strain
rep = args.rep 
hpi = args.hpi
fov1_idx = args.fov1
fov2_idx = args.fov2
base_dir = "/data/influenza-genome-packaging/results/preprocessed"
path_fov1 = f"{base_dir}/{strain}/rep{rep}/{hpi}hpi/primary-fov_00{fov1_idx}.nc" if fov1_idx < 10 else f"{base_dir}/{strain}/rep{rep}/{hpi}hpi/primary-fov_0{fov1_idx}.nc"
path_fov2 = f"{base_dir}/{strain}/rep{rep}/{hpi}hpi/primary-fov_00{fov2_idx}.nc" if fov2_idx < 10 else f"{base_dir}/{strain}/rep{rep}/{hpi}hpi/primary-fov_0{fov2_idx}.nc"

## Functions
def get_range(img):
    x_range = img.xc.values[[0, -1]]
    y_range = img.yc.values[[0, -1]]
    
    return x_range, y_range

def overlap_1d(range1, range2):
    overlap = [max(range1[0], range2[0]), min(range1[1], range2[1])]
    if(overlap[1]>overlap[0]):
        return overlap
    return None

def overlap_2d(ranges1, ranges2):
    xrange1, yrange1 = ranges1
    xrange2, yrange2 = ranges2
    return overlap_1d(xrange1, xrange2), overlap_1d(yrange1, yrange2)

## Load images and MIP
fov1 = xr.open_dataset(path_fov1).__xarray_dataarray_variable__
fov2 = xr.open_dataset(path_fov2).__xarray_dataarray_variable__
mp1 = fov1.max(dim=["r", "c"]).squeeze()
mp2 = fov2.max(dim=["r", "c"]).squeeze()

## Calculate pixel size
px_size_x1 = fov1.xc.values[1] - fov1.xc.values[0]
px_size_y1 = fov1.yc.values[1] - fov1.yc.values[0]
px_size_x2 = fov2.xc.values[1] - fov2.xc.values[0]
px_size_y2 = fov2.yc.values[1] - fov2.yc.values[0]

if px_size_x1 != px_size_x2 or px_size_y1 != px_size_y2 or px_size_x1 != px_size_x2:
    print("WARNING: pixel sizes differ between fov1 and fov2")

## Add channel and round images
viewer = napari.Viewer()

translates = [[0, fov1.yc[0], fov1.xc[0]], [0, fov2.yc[0], fov2.xc[0]]]

for fov, fov_idx in [(fov1, 1), (fov2, 2)]:
    scale = [1, px_size_y1, px_size_x1]

    for i,b,c in zip(range(3,-1,-1), "AGTC"[::-1], ["magenta", "green", "yellow", "red"][::-1]):
        viewer.add_image(fov[:,i], name=b + f" fov{fov_idx}", colormap=c, blending="additive", contrast_limits=[0, 0.4], visible=False, translate=translates[fov_idx-1], scale=scale)

## Add MIPs of FOV1 and FOV2
for mp, mp_idx, c in [(mp1, 1, "green"), (mp2, 2, "red")]:
    x0 = mp.xc[0]
    y0 = mp.yc[0]

    translate = [y0, x0]
    scale = [px_size_y1, px_size_x1]
    viewer.add_image(mp, name=f" mp{mp_idx}", colormap=c, blending="additive", contrast_limits=[0, 0.4], visible=True, translate=translate, scale=scale)

## Add Overlap Rectangle
# Calculate FOV ranges, widths and heights
range_fov1, range_fov2 = get_range(fov1), get_range(fov2)
width_fov1 = max(range_fov1[0]) - min(range_fov1[0])
height_fov1 = max(range_fov1[1]) - min(range_fov1[1])
width_fov2 = max(range_fov2[0]) - min(range_fov2[0])
height_fov2 = max(range_fov2[1]) - min(range_fov2[1])

# Calculate overlap range, width and height
overlap = overlap_2d(range_fov1, range_fov2)
x_overlap, y_overlap = overlap

width_overlap = max(x_overlap) - min(x_overlap)
height_overlap = max(y_overlap) - min(y_overlap)
width_overlap_px = width_overlap / px_size_x1
height_overlap_px = height_overlap / px_size_y1

# Add Rectangle
rect_ov = [[
    (0, 0),
    (height_overlap_px, width_overlap_px)
]]

viewer.add_shapes(
    rect_ov,
    shape_type="rectangle",
    edge_color="yellow",
    edge_width=5,
    face_color="transparent",
    opacity=0.5,
    name="Overlap Area",
    scale=[px_size_y1, px_size_x1], 
    translate=[fov1.yc.values[0] + height_fov1 - height_overlap, fov1.xc.values[0] + width_fov1 - width_overlap]
)

## Registration
# Cut out the overlap area of both MIPs
for i, fov, fov_idx, c in [(1, mp1, fov1_idx, "green"), (2, mp2, fov2_idx, "red")]:
    x_idx = np.where((fov.xc.values >= x_overlap[0]) & (fov.xc.values <= x_overlap[1]))[0]
    y_idx = np.where((fov.yc.values >= y_overlap[0]) & (fov.yc.values <= y_overlap[1]))[0]

    fov_ov = fov.sel(
        x = slice(x_idx[0], x_idx[-1]),
        y = slice(y_idx[0], y_idx[-1])
    )

    viewer.add_image(fov_ov, translate=[fov_ov.yc.values[0], fov_ov.xc.values[0]], scale=scale, name=f"FOV{i} overlap", colormap=c, blending="additive")



napari.run()