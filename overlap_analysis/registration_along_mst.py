from glob import glob
import xarray as xr
from itertools import combinations
import json
import pandas as pd
import os
import networkx as nx
from functions import common_roi, imagestack_from_netcdf, get_ranges_from_netcdf, get_range, overlap_1d, overlap_2d, overlap_area, overlap_length, crop_overlap, registration

strains = ["PR8", "Nepal", "PR8_Nepal"]#, "mutants/HA", "mutants/M1", "mutants/M3", "mutants/M4", "mutants/M5", 
           #"mutants/NP", "mutants/PA1", "mutants/PA2", "mutants/PB1", "mutants/PB2", "mutants/PSI2a", "mutants/PSIall", "mutants/WT"]
regs = []

for strain in strains:
    # if strain != "PR8":
    #     continue
    prefix = f"/data/influenza-genome-packaging/results/preprocessed/{strain}"
    reps = len(os.listdir(prefix))

    for rep in range(reps):
        # if rep != 2:
        #     continue
        hpis = sorted([int(hpi[0]) for hpi in os.listdir(f"{prefix}/rep{rep}")])

        for hpi in hpis:
            # if hpi != 1:
            #     continue
            base_dir = f"/data/influenza-genome-packaging/results/preprocessed/{strain}/rep{rep}/{hpi}hpi"
            fovs = sorted(glob(f"{base_dir}/primary-fov_*.nc"))

            # calculate fov ranges
            ranges = [get_ranges_from_netcdf(f) for f in fovs]

            # calculate Maximum Spanning Tree (MST)
            G = nx.Graph()

            for i1,i2 in combinations(range(len(ranges)), 2):
                r1 = ranges[i1]
                r2 = ranges[i2]
                area = overlap_area(overlap_2d(r1,r2))
                if area>0:
                    G.add_edge(i1, i2, weight=area)

            mst = nx.maximum_spanning_tree(G)
            edges = list(mst.edges)

            # Registration of all edges
            for source, target in edges:
                print(f"strain: {strain}, rep: {rep}, hpi: {hpi}, source: {source}, target: {target}")
                source_path = f"{base_dir}/primary-fov_0{source}.nc" if source >= 10 else f"{base_dir}/primary-fov_00{source}.nc"
                target_path = f"{base_dir}/primary-fov_0{target}.nc" if target >= 10 else f"{base_dir}/primary-fov_00{target}.nc"

                with xr.open_dataset(source_path, engine="netcdf4") as ds1, xr.open_dataset(target_path, engine="netcdf4") as ds2:
                    source_img = ds1.__xarray_dataarray_variable__.load()
                    target_img = ds2.__xarray_dataarray_variable__.load()

                source_img = xr.open_dataset(source_path).__xarray_dataarray_variable__
                target_img = xr.open_dataset(target_path).__xarray_dataarray_variable__
                source_max = source_img.max(dim=["r", "c"]).squeeze()
                target_max = target_img.max(dim=["r", "c"]).squeeze()

                px_size_x1 = source_img.xc.values[1] - source_img.xc.values[0]
                px_size_y1 = source_img.yc.values[1] - source_img.yc.values[0]
                px_size_x2 = target_img.xc.values[1] - target_img.xc.values[0]
                px_size_y2 = target_img.yc.values[1] - target_img.yc.values[0]

                warning = None
                if px_size_x1 != px_size_x2 or px_size_y1 != px_size_y2 or px_size_x1 != px_size_x2:
                    warning = "pixel sizes differ between source_img and target_img"


                range_source, range_target = get_range(source_img), get_range(target_img)
                mp_ov1, mp_ov2 = crop_overlap(source_max, target_max, range_source, range_target)

                # Warnings
                if mp_ov1.size == 0 or mp_ov2.size == 0:
                    warning = "overlap is smaller than pixel size"  
                    regs.append({"strain": strain, "rep": rep, "hpi": hpi, "source": source, "target":target, "tx": None, "ty": None, "warning": warning})
                    continue 
                
                if len(mp_ov1.xc.values) <= 1 or len(mp_ov2.xc.values) <= 1:
                    warning = "x-overlap is just one pixel"
                    regs.append({"strain": strain, "rep": rep, "hpi": hpi, "source": source, "target":target, "tx": None, "ty": None, "warning": warning})
                    continue

                if len(mp_ov1.yc.values) <= 1 or len(mp_ov2.yc.values) <= 1:
                    warning = "y-overlap is just one pixel"
                    regs.append({"strain": strain, "rep": rep, "hpi": hpi, "source": source, "target":target, "tx": None, "ty": None, "warning": warning})
                    continue

                tx_px, ty_px = registration(mp_ov1, mp_ov2, "registration_py.json")
                tx = tx_px * px_size_x1
                ty = ty_px * px_size_y1
                
                regs.append({"strain": strain, "rep": rep, "hpi": hpi, "source": source, "target":target, "tx": tx, "ty": ty, "warning": warning})

df = pd.DataFrame(regs)
df.to_csv("edges_and_registration.csv", index=False)
os.remove("registration_py.json")

