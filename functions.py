import pandas as pd
import os
import itertools
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from scipy.spatial import distance
import re
from starfish import Codebook
from starfish.core.spots.DetectPixels.pixel_spot_decoder import MultiBarcodePixelSpotDecoder
from starfish.core.imagestack.imagestack import ImageStack
from starfish.types import Axes
from scipy.spatial.distance import cdist
from scipy.stats import norm, ttest_1samp

total_reps = 5
total_hpis = 9
total_fovs = 10

barcodes_3nt = ["AGC","CGC","GGC","TAC","TCC","TTC","TGA","TGT"]
barcode_dict = {"PB2": "AGC", "PB1": "CGC", "PA": "GGC", "HA": "TAC", "NP": "TCC", "NA": "TTC", "M" :"TGA", "NS": "TGT"}
filenames = os.listdir("manual_annotation_decoded/")

def get_rep_hpi_fov(filename:str) -> tuple[int, int, int]:
    """This function extracts the replicate, hour post infection (hpi) and field of view (fov) from the filename.

    Args:
        filename (str): filename

    Returns:
        tuple[int, int, int]: Replicate, hour post infection, field of view
    """
    file_segments = filename.split("_")
    rep = file_segments[0][-1]
    hpi = file_segments[1][-1]
    fov = file_segments[2][-5]
    return int(rep), int(hpi), int(fov)

def get_bases_from_segments(segments:str) -> str:
    """Returns the bases of a viral segment.

    Args:
        segments (str): viral segment

    Returns:
        str: bases
    """
    if segments == "None":
        return 
    if "invalid" in segments or "missing" in segments:
        return "invalid"
    segments = segments.split(",")
    bases_man = ""
    bases_pipe = ""
    bases3 = ""
    for seg in segments:
        if str(barcode_dict[seg][0]) not in bases_man:
            bases_man += barcode_dict[seg][0]
        if str(barcode_dict[seg][1]) not in bases_pipe:
            bases_pipe += barcode_dict[seg][1]
        if str(barcode_dict[seg][2]) not in bases3:
            bases3 += barcode_dict[seg][2]
    return bases_man + "," + bases_pipe + "," + bases3

def get_spots(df_man, df_pipe, rep:int, hpi:int, fov:int):
    """Extracting specific spots for rep, hpi and fov out of manual annotated spots and pipeline annotated spots dataframes.

    Args:
        df_man (pd.Dataframe): Dataframe manual annotated spots
        df_pipe (pd.Dataframe): Dataframe pipeline annotated spots
        rep (int): Replication
        hpi (int): hour post infection
        fov (int): field of view

    Returns:
        pd.Dataframe: Dataframes with extracted spots for manual annotation and pipeline annotation
    """
    manual_spots_file = df_man[(df_man["rep"] == rep) & (df_man["hpi"] == hpi) & (df_man["fov"] == fov)]
    pipeline_spots_file = df_pipe[(df_pipe["rep"] == rep) & (df_pipe["hpi"] == hpi) & (df_pipe["fov"] == fov)]
    return manual_spots_file, pipeline_spots_file

def get_key_by_value(dictionary:dict, search:str) -> str:
    """Value based search for key in a dictionary.

    Args:
        dictionary (dict): dictionary
        search (str): value

    Returns:
        str: key
    """
    for key, value in dictionary.items():
        if value == search:
            return key
    return None 


def match_spots(list1:list, list2:list, max_distance:int) -> list:
    """Match spots from two lists based on coordinates. 

    Args:
        list1 (list): list of spots 1
        list2 (list): list of spots 2
        max_distance (int): maximal euclidean distance for matching

    Returns:
        list: list of matched spots
    """
    matched_spots = []
    used_spots1 = set()
    used_spots2 = set()

    for i, spot1 in enumerate(list1):
        closest_spot = None
        min_dist = float('inf')

        for j, spot2 in enumerate(list2):
            if j in used_spots2:
                continue

            dist = distance.euclidean(spot1, spot2)
            if dist < min_dist and dist < max_distance:
                min_dist = dist
                closest_spot = j

        if closest_spot is not None:
            matched_spots.append((i, closest_spot))
            used_spots1.add(i)
            used_spots2.add(closest_spot)

    return matched_spots

def get_unique_bases(bases_man:str, bases_pipe:str):
    """Takes two strings of bases as inputs and compares the differences of the bases of each run

    Args:
        bases_man (str): Bases of manual annotation as comma-seperated string
        bases_pipe (str): Bases of pipeline annotation as comma-seperated string

    Returns:
        list: Two lists of lists with the unique bases of each run, one for the manual annotated bases and one for the pipeline annotated bases
    """
    bases_man_split = bases_man.split(",")
    bases_pipe_split = bases_pipe.split(",")
    unique_bases_man = []
    unique_bases_pipe = []

    for run in range(3):
        bases_man_run = list(bases_man_split[run])
        bases_pipe_run = list(bases_pipe_split[run])
        unique_bases_man.append(list(set(bases_man_run) - set(bases_pipe_run)))
        unique_bases_pipe.append(list(set(bases_pipe_run) - set(bases_man_run)))
    
    return unique_bases_man, unique_bases_pipe

def common_roi(transforms_list):
    # crop the images to the common region
    translations = np.array([x[2].translation for x in transforms_list.transforms])
    lower = -np.floor(translations.min(axis=0)).astype(int)
    upper = -np.ceil(translations.max(axis=0)).astype(int)
    # the `or None` is to handle the case where no cropping is necessary at the upper bound (see: https://stackoverflow.com/a/11337953/4969760)
    return slice(lower[0],upper[0] or None), slice(lower[1],upper[1] or None)#Liste des Prozentualen Verlustes von Spots in der min_sigma range 1.0 bis 2.0:

def count_words(list: list, words_to_count: list) -> dict:
    """Counting words in an list

    Args:
        list (list): List with words
        words_to_count (list): List of words that should be counted

    Returns:
        dict: Dictionary with the word an the count
    """
    counts = {}
    for word in list:
        words = word.split(",")
        
        for w in words:
            if w in words_to_count:
                if f"{w}" in counts.keys():
                    counts[f"{w}"] += 1
                else:
                    counts[f"{w}"] = 1

    return counts

def iterative_pbd(image, start_dist_thresh=0.25, dist_thresh_interval=0.1, start_mag_thresh=0.30, mag_thresh_interval=0.15, intervals=2):
    single_codebook = Codebook.open_json('codebook_single_segment.json')
    decoded_pixels_all = pd.DataFrame()

    cumulative_mask = np.zeros(image.xarray.shape, dtype=bool)
    current_image_stack = image

    for iter in range(intervals):
        ## Decoding
        dist_thresh = start_dist_thresh + dist_thresh_interval * iter
        mag_thresh = start_mag_thresh - mag_thresh_interval * iter

        decoder = MultiBarcodePixelSpotDecoder(
            codebook=single_codebook,
            metric="cosine",
            norm_order=2,
            distance_threshold=dist_thresh,
            magnitude_threshold=mag_thresh,
        )

        decoded_pixels = decoder.run(current_image_stack)
        decoded_pixels_df = decoded_pixels.to_features_dataframe()
        decoded_pixels_df["iter"] = iter
        decoded_pixels_all = pd.concat([decoded_pixels_all, decoded_pixels_df])

        # Mask
        iter_mask = np.zeros(image.xarray.shape, dtype=bool)

        for _, row in decoded_pixels_df.iterrows():
            x, y = int(row["x"]), int(row["y"])
            r =  1
            y_min, y_max = max(0, y - r), min(iter_mask.shape[-2], y + r + 1)
            x_min, x_max = max(0, x - r), min(iter_mask.shape[-1], x + r + 1)
            iter_mask[..., y_min:y_max, x_min:x_max] = True

        cumulative_mask |= iter_mask
        masked_array = image.xarray.values.copy()
        masked_array[cumulative_mask] = 0

        current_image_stack = ImageStack.from_numpy(
            masked_array,
            index_labels={
                Axes.ROUND: list(range(masked_array.shape[0])),
                Axes.CH: list(range(masked_array.shape[1])),
                Axes.ZPLANE: list(range(masked_array.shape[2])),
            }
        )

    return decoded_pixels_all

def get_distance_from_bases_overlap(bases_man, bases_pipe):
    bases_man = str(bases_man)
    bases_pipe = str(bases_pipe)
    distance = 0

    if (bases_man == "invalid" and bases_pipe == "None") or (bases_man == "None" and bases_pipe == "invalid"):  ## None und invalid
        return "None/invalid"
    
    if bases_man == "invalid" and bases_pipe == "invalid":
        return "both invalid"
    
    if (bases_man == "invalid" and bases_pipe != "invalid") or (bases_man != "invalid" and bases_pipe == "invalid"): ## Invalid und valid
        return "invalid/valid"

    if (bases_man == "None" and bases_pipe != "None") or (bases_man != "None" and bases_pipe == "None"): ## None und valid
        bases = bases_man.split(",") if bases_man != "None" else bases_pipe.split(",")

        for bases_run in bases:
            distance += len(bases_run)
        return "None/valid"
    
    if bases_man not in ["None", "invalid"] and bases_pipe not in ["None", "invalid"]:
        bases_man = bases_man.split(",")
        bases_pipe = bases_pipe.split(",")

        for i in range(len(bases_man)):
            all_bases_run = bases_man[i] + bases_pipe[i]
            counts = Counter(all_bases_run)
            unique_counts = sum(1 for char in counts if counts[char] == 1)
            distance += unique_counts
        return distance
    
    return "No case"
    
def get_distance_from_bases(bases_man:str, bases_pipe:str) -> str:
    """Calculate barcode distance as number of different bases.

    Args:
        bases_man (str): bases from manual annotation
        bases_pipe (str): bases from pipeline annotation

    Returns:
        str: barcode distance
    """
    bases_man = str(bases_man)
    bases_pipe = str(bases_pipe)
    distance = 0

    if (re.search(r"[01]", bases_man) and bases_pipe == "None") or (bases_man == "None" and bases_pipe == "invalid"):  ## None und invalid
        return "both invalid"
    
    if re.search(r"[01]", bases_man) and bases_pipe == "invalid": ## invalid invalid
        return "both invalid"
    
    if (re.search(r"[01]", bases_man) and bases_pipe != "invalid|None") or (bases_man != re.search(r"[01]", bases_man) and bases_pipe == "invalid"): ## Invalid und valid
        return "invalid/valid"
    
    if (bases_man == "None" and bases_pipe != "None") or (bases_man != "None" and bases_pipe == "None"): ## None und valid
        bases = bases_man.split(",") if bases_man != "None" else bases_pipe.split(",")
        for bases_run in bases:
            distance += len(bases_run)
        return "None/valid"
    
    if (bases_man != "None" and not re.search(r"[01]", bases_man)) and bases_pipe not in ["None", "invalid"]: ## valid valid
        bases_man = bases_man.split(",")
        bases_pipe = bases_pipe.split(",")

        for i in range(len(bases_man)):
            all_bases_run = bases_man[i] + bases_pipe[i]
            counts = Counter(all_bases_run)
            unique_counts = sum(1 for char in counts if counts[char] == 1)
            distance += unique_counts
        return distance
    
    return "No case"

def pivot_replicates(df, index_cols, value_cols, rep_col="rep"):
    """
    Generates pivot tables for multiple experimental repeats (replicates) and combines them into a single DataFrame.

    Parameters:
    ----------
    df : pd.DataFrame
        The original DataFrame containing replicate data.
    index_cols : list of str
        Columns to be used as the index for the pivot table (e.g., ["hpi", "complex_size", "complex"]).
    value_cols : list of str
        Columns whose values across replicates should be represented as new columns.
    rep_col : str (default: "rep")
        The column that identifies the different replicates.

    Returns:
    -------
    pd.DataFrame
        A merged pivot table with columns like: value_rep0, value_rep1, ...
    """

    pivoted_tables = []

    for value in value_cols:
        pivot = df.pivot(index=index_cols, columns=rep_col, values=value)
        pivot.columns = [f"{value}_{rep_col}{col}" for col in pivot.columns]

        mean = np.mean(pivot.values, axis=1)
        std = np.std(pivot.values, axis=1, ddof=1)  # ddof=1 = Stichproben-Standardabweichung
        pivot[f"{value}_mean"] = mean
        pivot[f"{value}_std"] = std

        pivoted_tables.append(pivot)

    result = pd.concat(pivoted_tables, axis=1).reset_index()
    return result

def simulated_error(n_simulations, complexes_df, result_folder):
    df = complexes_df.loc[:, ["n_sum", "relative_abundance_mean", "normalized_complex_probability_mean"]]
    
    mean_relative_abundances = complexes_df["relative_abundance_mean"]
    mean_expected_probabilities = complexes_df["normalized_complex_probability_mean"]
    n_complexes = sum(complexes_df["n_sum"])

    mean_errors = []
    simulation_data = []

    if n_complexes == 0:
        print("no complexes")
        return None, None, None, None, None
    

    # Simulationen durchführen
    for i in range(n_simulations):
        simulated_counts = np.random.multinomial(n_complexes, mean_expected_probabilities)
        simulated_relative_abundance = simulated_counts / n_complexes
        mean_error = np.mean(np.abs(simulated_relative_abundance - mean_expected_probabilities))
        mean_errors.append(mean_error)
        simulation_data.append(pd.Series(simulated_counts, index=complexes_df.index, name=f"simulation{i}_count"))
    
    # Simulationsdaten hinzufügen
    df = pd.concat([df] + simulation_data, axis=1)

    #Mittelwert und Standardabweichung berechnen
    mean_mean_error = np.mean(mean_errors)
    std_mean_error = np.std(mean_errors)
    observed_mean_error = np.mean((np.abs(mean_relative_abundances - mean_expected_probabilities)))

    # Z-Score und zweiseitiger p-Wert für Vergleichszwecke
    z_score = (observed_mean_error - mean_mean_error) / std_mean_error
    p_value = 2 * (1 - norm.cdf(np.abs(z_score)))
    p_value_corrected = p_value * len(complexes_df) 

    # Einseitiger t-Test für observed_mean_error > mean_mean_error
    p_value_one_sided = p_value / 2 if observed_mean_error > mean_mean_error else 1

    # Speichern der Simulationsergebnisse
    df.to_csv(f"{result_folder}.csv", index=False)

    # Plotten der Ergebnisse
    plt.hist(mean_errors, bins=30, alpha=0.7, color="blue")
    plt.axvline(x=mean_mean_error, color="red", linestyle="--", label="simulated mean error")
    plt.axvline(x=observed_mean_error, color="green", linestyle="--", label="observed mean error")
    plt.xlabel("mean error")
    plt.ylabel("frequency")
    plt.title("Verteilung der mittleren Fehler aus Simulationen")
    plt.legend()
    plt.savefig(result_folder + ".png")
    plt.close()

    return mean_mean_error, observed_mean_error, p_value, p_value_corrected, p_value_one_sided