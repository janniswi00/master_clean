import os
import glob
import numpy as np
import argparse
import pandas as pd
from scipy.ndimage import gaussian_filter
import xarray as xr
from scipy import ndimage

filenames = sorted(os.listdir("manual_annotation_decoded/"))
parser = argparse.ArgumentParser()
parser.add_argument('--exp', type=str, help='name of experiment', default="experiment")
args = parser.parse_args()
experiment_name = args.exp

def get_rep_hpi_fov(filename):
    file_segments = filename.split("_")
    rep = file_segments[0][-1]
    hpi = file_segments[1][-1]
    fov = file_segments[2][-5]
    return rep, hpi, fov

def common_roi(transforms_list):
    # crop the images to the common region
    translations = np.array([x[2].translation for x in transforms_list.transforms])
    lower = -np.floor(translations.min(axis=0)).astype(int)
    upper = -np.ceil(translations.max(axis=0)).astype(int)
    # the `or None` is to handle the case where no cropping is necessary at the upper bound (see: https://stackoverflow.com/a/11337953/4969760)
    return slice(lower[0],upper[0] or None), slice(lower[1],upper[1] or None)


parameters = {
    "masking_radius": 4,
    "blob_det_max_sigma": 5,
    "blob_det_thresh": 0.045,
    "pbd_dist_thresh": 0.25,
    "pbd_mag_thresh": 0.15,
    "density_measurement_max_int_thresh": 0.25,
    "gaussian_sigma": 5
}
total_hpis = 9
total_reps = 5
img_features = pd.read_csv("data/image_features_with_predictions.csv")

for rep in range(total_reps):
    for hpi in range(total_hpis):
        data_folder = f"pipeline_output/final_pipe/3nt/PR8/rep{rep}/{hpi}hpi/"
        total_fovs = len(glob.glob(f"{data_folder}*parameters.csv"))

        for fov_idx in range(total_fovs):
            if os.path.exists(f"pipeline_output/{experiment_name}/3nt/PR8/rep{rep}/{hpi}hpi/fov_{fov_idx}_density_map.npy"):
                print(f"{data_folder}{fov_idx}, done")
                continue

            primary_ds = xr.open_dataset(f"{data_folder}fov_{fov_idx}_primary.nc")
            primary = primary_ds[list(primary_ds.data_vars)[0]]
            summed_img = primary.sum(dim=("r", "c"))
            img_2d = summed_img.values.squeeze()
            density_map = gaussian_filter(img_2d, sigma=parameters["gaussian_sigma"])

            percentile = img_features[img_features["file"] == f"rep{rep}_hpi{hpi}_fov{fov_idx}"].predicted_percentile.values[0]

            if hpi <= 3:
                percentile = 100

            if percentile > 100:
                percentile = 100

            print(f"rep{rep}_hpi{hpi}_fov{fov_idx}", percentile)

            threshold_percentile = np.percentile(density_map, percentile)
            mask_percentile = (density_map > threshold_percentile).astype(int)

            labeled_mask, num_features = ndimage.label(mask_percentile)
            sizes = ndimage.sum(mask_percentile, labeled_mask, range(1, num_features + 1))

            min_size = 200
            mask_filtered = np.zeros_like(mask_percentile, dtype=bool)

            for i, size in enumerate(sizes):
                if size >= min_size:
                    mask_filtered[labeled_mask == (i + 1)] = True

            dir = f"pipeline_output/{experiment_name}/3nt/PR8/rep{rep}/{hpi}hpi"
            os.makedirs(dir, exist_ok=True)
            prefix = f"{dir}/fov_{fov_idx}"
            pd.DataFrame(list(parameters.items()), columns=["parameters", "value"]).to_csv(f"{prefix}_parameters.csv", index=False)
            np.save(f"{prefix}_density_map.npy", density_map)
            np.save(f"{prefix}_density_mask.npy", mask_percentile)
            np.save(f"{prefix}_density_mask_filtered.npy", mask_filtered)
