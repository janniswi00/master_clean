from glob import glob
import numpy as np
import xarray as xr
from starfish.core.imagestack.parser.numpy import NumpyData
from starfish.core.imagestack.imagestack import ImageStack
from itertools import combinations
import json
from starfish.types import Axes
import matplotlib.pyplot as plt
from starfish.image import ApplyTransform, LearnTransform, Filter
import json

def common_roi(transforms_list):
    # crop the images to the common region
    translations = np.array([x[2].translation for x in transforms_list.transforms])
    lower = -np.floor(translations.min(axis=0)).astype(int)
    upper = -np.ceil(translations.max(axis=0)).astype(int)
    # the `or None` is to handle the case where no cropping is necessary at the upper bound (see: https://stackoverflow.com/a/11337953/4969760)
    return slice(lower[0],upper[0] or None), slice(lower[1],upper[1] or None)

def imagestack_from_netcdf(filename):
    imgs = xr.open_dataset(filename).__xarray_dataarray_variable__
    # using imgs.indexes does not work because the values in the dict are instances of class Index
    index_labels = {x:list(y) for x,y in dict(imgs.indexes).items()}
    image_stack = ImageStack.from_tile_collection_data(NumpyData(imgs, index_labels, imgs.coords))
    return image_stack

def get_ranges_from_netcdf(filename):
    img = imagestack_from_netcdf(filename)
    xrange = img.xarray.xc.values[[0,-1]]
    yrange = img.xarray.yc.values[[0,-1]]
    return xrange, yrange

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

def overlap_length(overlap):
    if not overlap:
        return 0
    return overlap[1]-overlap[0]

def overlap_area(overlaps):
    return overlap_length(overlaps[0]) * overlap_length(overlaps[1])

def crop_overlap(max1, max2, range1, range2):
    overlap = overlap_2d(range1, range2)
    x_overlap, y_overlap = overlap
    mps_ov = []

    for max in [max1, max2]:
        x_idx = np.where((max.xc.values >= x_overlap[0]) & (max.xc.values <= x_overlap[1]))[0]
        y_idx = np.where((max.yc.values >= y_overlap[0]) & (max.yc.values <= y_overlap[1]))[0]

        max_ov = max.sel(
            x = slice(x_idx[0], x_idx[-1]),
            y = slice(y_idx[0], y_idx[-1])
        )

        mps_ov.append(max_ov)

    return mps_ov

def registration(mp_ov1, mp_ov2):
    stack_np = np.stack([mp_ov1.values, mp_ov2.values], axis=0)
    stack_np = stack_np[:, np.newaxis, np.newaxis, :, :]
    stack_both = ImageStack.from_numpy(stack_np)

    learn_translation = LearnTransform.Translation(reference_stack=stack_both.sel({Axes.ROUND: 0}), axes=Axes.ROUND, upsampling=1000) ## Registration
    warp = ApplyTransform.Warp()
    transform_list = learn_translation.run(stack_both)
    registered_stack = warp.run(stack_both, transforms_list=transform_list)
    transform_list.to_json(f"registration2.json")

    with open(f"registration2.json", "r", encoding="utf-8") as file_reg:
        dict_reg = json.load(file_reg)

    tx = dict_reg["transforms_list"][1][2][0][2] * -1
    ty = dict_reg["transforms_list"][1][2][1][2] * -1
    
    return tx, ty