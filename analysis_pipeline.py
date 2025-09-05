#!/usr/bin/env python
# coding: utf-8

# # Starfish-based analysis pipeline for 3nt data
import os
import numpy as np
import skimage.io
import scipy.ndimage as ndi
from skimage.filters import threshold_otsu
from starfish import Experiment, Codebook
from starfish.types import Axes, Levels
from starfish.image import ApplyTransform, LearnTransform, Filter
from starfish.spots import FindSpots, DecodeSpots
from starfish.morphology import Binarize, Merge, Segment
import starfish.morphology # Filter is already imported from starfish.image
from starfish.core.spots.DetectPixels.pixel_spot_decoder import MultiBarcodePixelSpotDecoder
from tqdm import tqdm
import argparse
from sklearn.cluster import DBSCAN
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy import ndimage
from functions import iterative_pbd

parser = argparse.ArgumentParser()
parser.add_argument('--exp', type=str, help='name of the experiment', default="experiment")
args = parser.parse_args()
experiment_name = args.exp
img_features = pd.read_csv("image_features_with_predictions.csv")

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
    "gaussian_sigma": 5,
    "mask_filter_min_size": 200,
    "pbd_start_dist_thresh": 0.25,
    "pbd_dist_thresh_interval": 0.05, 
    "pbd_start_mag_thresh": 0.30, 
    "pbd_mag_thresh_interval": 0.075, 
    "pbd_intervals": 3,
}

total_reps = 5
total_hpis = 9

for rep in range(total_reps):
    for hpi in range(total_hpis):
        datafolder = f"/tank/s391697/in-situ-seq-influenza/data/starfish/3nt/PR8/rep{rep}/{hpi}hpi"
        exp = Experiment.from_json(os.path.join(datafolder, "experiment.json"))

        for fov_index, fov in tqdm(enumerate(exp.fovs())):
            if os.path.exists(f"pipeline_output/{experiment_name}/3nt/PR8/rep{rep}/{hpi}hpi/fov_{fov_index}_primary.nc"):
                print(f"rep{rep}/{hpi}hpi/fov_{fov_index} done")
                continue

            imgs = fov.get_image("primary")
            nuclei = fov.get_image("nuclei")
            bf = fov.get_image("brightfield")

            # Spot based registration
            spots = imgs.reduce({Axes.CH}, func="max")
            learn_translation = LearnTransform.Translation(reference_stack=spots.sel({Axes.ROUND: 0}), axes=Axes.ROUND, upsampling=1000)
            warp = ApplyTransform.Warp()
            transforms_list = learn_translation.run(spots)
            registered_imgs = warp.run(imgs, transforms_list=transforms_list)
            transforms_list = learn_translation.run(spots)
            registered_nuclei = warp.run(nuclei, transforms_list=transforms_list)
            transforms_list = learn_translation.run(spots)
            registered_bf = warp.run(bf, transforms_list=transforms_list)
            # re-create for saving as json
            transforms_list = learn_translation.run(spots)


            # ### Crop to common region
            xroi, yroi = common_roi(transforms_list)
            registered_imgs = registered_imgs.sel({Axes.X: xroi, Axes.Y: yroi})
            registered_nuclei = registered_nuclei.sel({Axes.X: xroi, Axes.Y: yroi})
            registered_bf = registered_bf.sel({Axes.X: xroi, Axes.Y: yroi})

            # ## Background Filtering
            masking_radius = parameters["masking_radius"]
            filt = Filter.WhiteTophat(masking_radius, is_volume=False)
            filtered = filt.run(registered_imgs, verbose=True, in_place=False)

            # ## Normalization
            cptz= Filter.ClipPercentileToZero(p_min=80, p_max=100, group_by={Axes.CH}, level_method=Levels.SCALE_BY_CHUNK)
            clipped_both_scaled = cptz.run(filtered, in_place=False)

            ## Density measurement
            primary = clipped_both_scaled.xarray
            summed_img = primary.sum(dim=("r", "c"))
            img_2d = summed_img.values.squeeze()
            density_map = gaussian_filter(img_2d, sigma=parameters["gaussian_sigma"])

            percentile = img_features[img_features["file"] == f"rep{rep}_hpi{hpi}_fov{fov_index}"].predicted_percentile.values[0]

            if hpi <= 3:
                percentile = 100

            if percentile > 100:
                percentile = 100

            threshold_percentile = np.percentile(density_map, percentile)
            mask_percentile = (density_map > threshold_percentile).astype(int)

            labeled_mask, num_features = ndimage.label(mask_percentile)
            sizes = ndimage.sum(mask_percentile, labeled_mask, range(1, num_features + 1))

            mask_filtered = np.zeros_like(mask_percentile, dtype=bool)

            for i, size in enumerate(sizes):
                if size >= parameters["mask_filter_min_size"]:
                    mask_filtered[labeled_mask == (i + 1)] = True

            # ## Spot Detection
            dots = clipped_both_scaled.reduce([Axes.CH, Axes.ROUND], "max")
            bd = FindSpots.BlobDetector(
                min_sigma=1,
                max_sigma=parameters["blob_det_max_sigma"],
                num_sigma=20,
                threshold=parameters["blob_det_thresh"],
                is_volume=True,
                measurement_type='max',
            )
            spots = bd.run(clipped_both_scaled, reference_image=dots)
            
            # ## Decode Spots spot-based
            single_codebook = Codebook.open_json('codebook_single_segment.json')

            multi_barcode_decoder = DecodeSpots.MultiBarcodeDecoder(
                codebook=single_codebook,
                max_distance=.5,
                min_intensity=.1,
                return_original_intensities=True
            )
            mbd_decoded_spots = multi_barcode_decoder.run(spots=spots)

            try:
                df_pixel = iterative_pbd(clipped_both_scaled, parameters["pbd_start_dist_thresh"],parameters["pbd_dist_thresh_interval"], parameters["pbd_start_mag_thresh"], parameters["pbd_mag_thresh_interval"], parameters["pbd_intervals"])
            except IndexError:
                  continue


            ## Clustering Pixel-based decoding
            eps_dist = 1
            min_pixels = 3

            coords = df_pixel[["x", "y"]].to_numpy()
            clustering = DBSCAN(eps=eps_dist, min_samples=min_pixels).fit(coords)
            df_pixel["spot_id"] = clustering.labels_
            df_pixel = df_pixel[df_pixel["spot_id"] != -1]

            spot_coords = df_pixel.groupby("spot_id")[["x", "y"]].mean().reset_index()
            df_pixel = df_pixel.drop(columns=[col for col in df_pixel.columns if col in ["x_center", "y_center"]])
            df_pixel = df_pixel.merge(spot_coords, on="spot_id", suffixes=("", "_center"))
            df_pixel["distance_to_center"] = np.sqrt((df_pixel['x'] - df_pixel['x_center'])**2 + (df_pixel['y'] - df_pixel['y_center'])**2)
            spot_radii = df_pixel.groupby("spot_id")["distance_to_center"].max().reset_index(name="radius")
            spot_coords = spot_coords.merge(spot_radii, on="spot_id")

            counts = df_pixel.groupby(["spot_id", "target"]).size().reset_index(name="count")
            total_counts = df_pixel.groupby("spot_id").size().reset_index(name="total")
            merged = counts.merge(total_counts, on="spot_id")
            merged["fraction"] = merged["count"] / merged["total"]

            filtered = merged[merged["fraction"] >= 0.1]

            spot_targets = (
                filtered[["spot_id", "target"]]
                .groupby("spot_id")["target"]
                .apply(lambda x: ",".join(sorted(set(x))))
                .reset_index()
            )
            feature_df_pbd = spot_coords.merge(spot_targets[["spot_id", "target"]], on="spot_id")

            feature_df_pbd["target"] = feature_df_pbd["target"].fillna("unassigned")

            spot_iters = df_pixel.groupby("spot_id")["iter"].unique().reset_index()
            spot_iters["iter"] = spot_iters["iter"].apply(lambda x: ",".join(map(str, sorted(set(x)))))

            feature_df_pbd = feature_df_pbd.merge(spot_iters, on="spot_id")
            
            # ## Segmentation

            # don't use starfish.image.Segment.Watershed because it uses hard coded thresholds for nucleus size (10,10000) which are not applicable here
            # see: https://spacetx-starfish.readthedocs.io/en/latest/gallery/tutorials/watershed_segmentation.html

            nuclei_red = registered_nuclei.reduce(dims=[Axes.ROUND, Axes.CH, Axes.ZPLANE], func="max")
            smoothed_nuclei = Filter.GaussianLowPass(sigma=5, is_volume=False).run(nuclei_red)
            dapi_thresh = threshold_otsu(smoothed_nuclei.xarray.values) # e.g. .048, binary mask for cell (nuclear) locations
            stain_thresh = 0.004  # binary mask for overall cells // binarization of stain
            min_dist = 56
            min_allowed_size = 1000
            max_allowed_size = 100000

            binarized_nuclei = Binarize.ThresholdBinarize(dapi_thresh).run(smoothed_nuclei)
            labeled_masks = starfish.morphology.Filter.MinDistanceLabel(min_dist, 1).run(binarized_nuclei)
            watershed_markers = starfish.morphology.Filter.AreaFilter(min_area=min_allowed_size, max_area=max_allowed_size).run(labeled_masks)

            thresholded_stain = Binarize.ThresholdBinarize(stain_thresh).run(nuclei_red)
            markers_and_stain = Merge.SimpleMerge().run([thresholded_stain, watershed_markers])
            watershed_mask = starfish.morphology.Filter.Reduce(
                "logical_or",
                lambda shape: np.zeros(shape=shape, dtype=bool)
            ).run(markers_and_stain)

            segmenter = Segment.WatershedSegment(connectivity=np.ones((1, 3, 3), dtype=bool))

            # masks is BinaryMaskCollection for downstream steps
            masks = segmenter.run(
                nuclei_red,
                watershed_markers,
                watershed_mask,
            )

            # ## Save Results
            feature_df_sbd = mbd_decoded_spots.to_features_dataframe()
            nucleus_mask = watershed_markers.to_label_image().xarray.squeeze(Axes.ZPLANE.value).values
            feature_df_sbd["nucleus"] = nucleus_mask[feature_df_sbd.y, feature_df_sbd.x]>0
            distance_to_nucleus = ndi.distance_transform_edt(~(nucleus_mask>0))-ndi.distance_transform_edt(nucleus_mask>0)
            feature_df_sbd["nucleus_dist"] = distance_to_nucleus[feature_df_sbd.y, feature_df_sbd.x]
            cell_mask = masks.to_label_image().xarray.squeeze(Axes.ZPLANE.value).values
            feature_df_sbd["cell"] = cell_mask[feature_df_sbd.y, feature_df_sbd.x]
            dist_to_boundary = ndi.distance_transform_edt(~skimage.segmentation.find_boundaries(cell_mask))
            feature_df_sbd["boundary_dist"] = dist_to_boundary[feature_df_sbd.y, feature_df_sbd.x]
            border_cells = np.unique(np.concatenate([cell_mask[:,0], cell_mask[:,-1], cell_mask[0,:], cell_mask[-1,:]]))
            feature_df_sbd["border_cell"] = [c in border_cells for c in feature_df_sbd["cell"]]
            feature_df_sbd["densed_region"] = mask_percentile[feature_df_sbd.y, feature_df_sbd.x].astype(bool)
            feature_df_sbd["densed_region_filtered"] = mask_filtered[feature_df_sbd.y, feature_df_sbd.x]

            feature_df_pbd["x"] = np.round(feature_df_pbd["x"], 0).astype(int)
            feature_df_pbd["y"] = np.round(feature_df_pbd["y"], 0).astype(int)
            feature_df_pbd["nucleus"] = nucleus_mask[feature_df_pbd.y, feature_df_pbd.x]>0
            feature_df_pbd["nucleus_dist"] = distance_to_nucleus[feature_df_pbd.y, feature_df_pbd.x]
            feature_df_pbd["cell"] = cell_mask[feature_df_pbd.y, feature_df_pbd.x]
            feature_df_pbd["boundary_dist"] = dist_to_boundary[feature_df_pbd.y, feature_df_pbd.x]
            feature_df_pbd["border_cell"] = [c in border_cells for c in feature_df_pbd["cell"]]
            feature_df_pbd["densed_region"] = mask_percentile[feature_df_pbd.y, feature_df_pbd.x].astype(bool)
            feature_df_pbd["densed_region_filtered"] = mask_filtered[feature_df_pbd.y, feature_df_pbd.x]

            dir = f"pipeline_output/{experiment_name}/3nt/PR8/rep{rep}/{hpi}hpi"
            os.makedirs(dir, exist_ok=True)
            prefix = f"{dir}/fov_{fov_index}"
            pd.DataFrame(list(parameters.items()), columns=["parameters", "value"]).to_csv(f"{prefix}_parameters.csv", index=False)
            skimage.io.imsave(f"{prefix}_nuclei.png", nucleus_mask)
            skimage.io.imsave(f"{prefix}_cells.png", cell_mask)
            feature_df_sbd.to_csv(f"{prefix}_spots_sbd.csv", na_rep=None, index=False)
            feature_df_pbd.to_csv(f"{prefix}_spots_pbd.csv", index=False)
            mbd_decoded_spots.to_netcdf(f"{prefix}_spots.nc")
            clipped_both_scaled.xarray.to_netcdf(f"{prefix}_primary.nc")
            transforms_list.to_json(f"{prefix}_registration.json")

            np.save(f"{prefix}_density_map.npy", density_map)
            np.save(f"{prefix}_density_mask.npy", mask_percentile)
            np.save(f"{prefix}_density_mask_filtered.npy", mask_filtered)
            print(f"done rep: {rep}, hpi: {hpi}, fov: {fov_index}")
