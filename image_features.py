import numpy as np
import pandas as pd
import glob
import xarray as xr 
from scipy.stats import skew, kurtosis, entropy

img_features = pd.DataFrame()
total_hpis = 9
total_reps = 5
percentile_annot = pd.read_csv("data/manual_percentiles.csv", sep=";")

for rep in range(total_reps):
    for hpi in range(total_hpis):
        data_folder = f"final_pipe/3nt/PR8/rep{rep}/{hpi}hpi/"
        total_fovs = len(glob.glob(f"{data_folder}*parameters.csv"))

        for fov_idx in range(total_fovs):

            ## load image
            primary_ds = xr.open_dataset(f"{data_folder}fov_{fov_idx}_primary.nc")
            primary = primary_ds[list(primary_ds.data_vars)[0]]
            summed_img = primary.sum(dim=("r", "c"))
            img_2d = summed_img.values.squeeze()
            manual_percentile_row = percentile_annot[(percentile_annot["rep"] == rep) & (percentile_annot["hpi"] == hpi) & (percentile_annot["fov"] == fov_idx)]
            if len(manual_percentile_row) > 0:
                manual_percentile = manual_percentile_row["percentile"].values[0] 
            else:
                manual_percentile = None

            ## calculate features
            hist, _ = np.histogram(img_2d, bins=20, range = (img_2d.min(), img_2d.max()), density=True)
            hist_features = {f"hist_bin_{i}": [v] for i, v in enumerate(hist)}
            p25, p50, p75, p90 = np.percentile(img_2d, [25, 50, 75, 90])
            bottom_90 = img_2d[img_2d <= np.percentile(img_2d, 90)]
            top_10 = img_2d[img_2d > np.percentile(img_2d, 90)]
            bottom_mean = bottom_90.mean() if bottom_90.size > 0 else 1e-6
            top_mean = top_10.mean() if top_10.size > 0 else 1e-6
            
            df = pd.DataFrame({
                "file": [f"rep{rep}_hpi{hpi}_fov{fov_idx}"],
                "min_intensity": [np.min(img_2d)],
                "max_intensity": [np.max(img_2d)],
                "mean_intensity": [np.mean(img_2d)],
                "std_intensity": [np.std(img_2d)],
                "skewness": [skew(img_2d.flatten())],
                "kurtosis": [kurtosis(img_2d.flatten())],
                "percentile_25": [p25],
                "percentile_50": [p50],
                "percentile_75": [p75],
                "percentile_90": [p90],
                "entropy": [entropy(hist + 1e-10)],
                "sbr": [top_mean / bottom_mean],
                "manual_percentile": [manual_percentile],
                **hist_features
            })

            img_features = pd.concat([img_features, df])

img_features.to_csv("data/image_features.csv", index=False)