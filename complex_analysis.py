import glob
import pandas as pd
import numpy as np
from scipy import stats

total_hpis = 9
total_reps = 5

paths = [f"/home/s391697/f2/complex_likelihood/results/complex_counts_all/complex_counts_rep{i}.xlsx" for i in range(total_reps)]
gene_names = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]

def calculate_complex_likelihood(paths_complex_counts, path_results, rank, nucleus, cytosol):
    df_all = []
    for rep, path in enumerate(paths_complex_counts):  ## für alle reps
        complex_counts_rep = pd.read_csv(path, keep_default_na=False)

        for hpi in range(total_hpis):  ## für alle hpis eines reps
            paths_fovs_csv = sorted(glob.glob(f"/home/s391697/f2/data/analysis/starfish/3nt/PR8/rep{rep}/{hpi}hpi/*spots.csv"))
            dfs_rep_hpi = []

            for fov_csv in paths_fovs_csv:  ## für alle fovs eines hpi eines reps
                df = pd.read_csv(fov_csv, keep_default_na=False)
                dfs_rep_hpi.append(df)
            
            combined_dfs_hpi = pd.concat(dfs_rep_hpi, ignore_index=True) ## alle spots eines hpis eines reps
            if nucleus == True:
                combined_dfs_hpi = combined_dfs_hpi[combined_dfs_hpi["nucleus"]==True]
            if cytosol == True:
                combined_dfs_hpi = combined_dfs_hpi[combined_dfs_hpi["nucleus"]==False]
            combined_dfs_hpi = combined_dfs_hpi[["target"]] ## nur die spalte mit gennamen
            combined_dfs_hpi = combined_dfs_hpi[~combined_dfs_hpi["target"].str.contains("invalid")]
            combined_dfs_hpi = combined_dfs_hpi[~combined_dfs_hpi["target"].str.contains("missing")]
    

            ## Berechnung der Häufigkeit jedes Segments einer hpi eines reps
            genes_dictionary = []

            for i in range(len(gene_names)):
                gene_abundance = (combined_dfs_hpi["target"].str.contains(gene_names[i])).sum()
                dictionary = {"gene": gene_names[i], "abundance": gene_abundance}
                genes_dictionary.append(dictionary)
            
            genes_dictionary_df = pd.DataFrame(genes_dictionary)

            ## Berechnung der relativen Häufigkeit jedes Segments einer hpi eines reps
            for i in range(len(genes_dictionary_df)):
                relative_abundance = genes_dictionary_df.iloc[i, 1] / genes_dictionary_df["abundance"].sum()
                genes_dictionary[i]["relative_abundance"] = relative_abundance
            
            gene_abundances = pd.DataFrame(genes_dictionary)


            ## Berechnung der Komplexwahrscheinlichkeit eines reps einer hpi
            complex_counts_rep_hpi = complex_counts_rep[complex_counts_rep["hpi"] == hpi].copy()
            complex_counts_rep_hpi.loc[:, "complex"] = complex_counts_rep_hpi["complex"].astype(str)
            if nucleus == True:
                complex_counts_rep_hpi = complex_counts_rep_hpi[complex_counts_rep_hpi["nucleus"] == True]
            if cytosol == True:
                complex_counts_rep_hpi = complex_counts_rep_hpi[complex_counts_rep_hpi["nucleus"] == False]

            ##relative Häufigkeit der Komplexe innerhalb der ranks
            complex_counts_rep_hpi["relative_abundance"] = (complex_counts_rep_hpi.groupby("complex_size")["n"].transform(lambda x: x / x.sum()))

            complex_probabilities = []
            for i in range(len(complex_counts_rep_hpi)):
                complex_genes = complex_counts_rep_hpi.iloc[i]["complex"].split(",")
                complex_probability = 1

                for gene in complex_genes:
                    relative_abundance = gene_abundances.loc[gene_abundances["gene"] == gene, "relative_abundance"].values
                    if relative_abundance.size > 0:
                        complex_probability *= relative_abundance[0]
                    else:
                        complex_probability = 0
                
                complex_probabilities.append(complex_probability)
            
            complex_counts_rep_hpi["complex probability"] = complex_probabilities

            complex_counts_rep_hpi_rank2 = complex_counts_rep_hpi[complex_counts_rep_hpi["complex_size"] == rank].copy() ## relative complex probability rank
            normalized_complex_probabilities = []
            sum_of_probabilities = complex_counts_rep_hpi_rank2["complex probability"].sum()

            for i in range(len(complex_counts_rep_hpi_rank2)):
                normalized_complex_probability = complex_counts_rep_hpi_rank2["complex probability"].iloc[i] / sum_of_probabilities
                normalized_complex_probabilities.append(normalized_complex_probability)
            
            complex_counts_rep_hpi_rank2["normalized_complex_probability"] = normalized_complex_probabilities
            complex_counts_rep_hpi_rank2["rep"] = rep  ## rep hinzufügen
            complex_counts_rep_hpi_rank2 = complex_counts_rep_hpi_rank2[["rep"] + [col for col in complex_counts_rep_hpi_rank2.columns if col != "rep"]]
            df_all.append(complex_counts_rep_hpi_rank2)

    final_df = pd.concat(df_all, ignore_index=True)
    final_df.to_excel(path_results + "complex_likelihood_all.csv", index=False)

    ## Create excel sheets with mean likelihoods across hpis
    for hpi in range(total_hpis):
        hpi_results = pd.DataFrame(columns=["hpi","complex","relative abundance rep0","relative abundance rep1","relative abundance rep2","relative abundance rep3","relative abundance rep4",
                                            "mean relative abundance", "sd relative abundance", "normalized complex probability rep0","normalized complex probability rep1",
                                            "normalized complex probability rep2","normalized complex probability rep3","normalized complex probability rep4",
                                            "mean normalized complex probability", "sd normalized complex probability", "ttest t_stat", "ttest p_value"])
        complexes_hpi_df = final_df[final_df["hpi"] == hpi]  ## filtern nach hpi
        complexes_hpi = complexes_hpi_df["complex"].unique()  ## einzigartige Komplexe ermitteln

        for complex in complexes_hpi:  ## berechnung von mittelwert und stabw
            complex_hpi = complexes_hpi_df[complexes_hpi_df["complex"] == complex]
            detected_in_reps = complex_hpi["rep"].unique()
            relative_abundance_values = []
            normalized_complex_probability_values = []
            complex_counts = []
            
            for rep in range(total_reps):  ## reps in denen kein Komplex detektiert wurde gleich 0 setzen
                if rep in detected_in_reps:
                    relative_abundance_values.append(complex_hpi[complex_hpi["rep"]==rep]["relative_abundance"].values[0])
                    normalized_complex_probability_values.append(complex_hpi[complex_hpi["rep"]==rep]["normalized_complex_probability"].values[0])
                    complex_counts.append(complex_hpi[complex_hpi["rep"]==rep]["n"].values[0])
                else:
                    relative_abundance_values.append(0)
                    normalized_complex_probability_values.append(0)
            

            t_stat, p_value = stats.ttest_ind(relative_abundance_values, normalized_complex_probability_values)

            complex_hpi_df = pd.DataFrame({
                "hpi": [hpi],
                "complex": [complex],
                "complex_count": [sum(complex_counts)],
                "relative abundance rep0": [relative_abundance_values[0]],
                "relative abundance rep1": [relative_abundance_values[1]],
                "relative abundance rep2": [relative_abundance_values[2]],
                "relative abundance rep3": [relative_abundance_values[3]],
                "relative abundance rep4": [relative_abundance_values[4]],
                "mean relative abundance": [np.mean(relative_abundance_values)],
                "sd relative abundance": [np.std(relative_abundance_values)],
                "normalized complex probability rep0": [normalized_complex_probability_values[0]],
                "normalized complex probability rep1": [normalized_complex_probability_values[1]],
                "normalized complex probability rep2": [normalized_complex_probability_values[2]],
                "normalized complex probability rep3": [normalized_complex_probability_values[3]],
                "normalized complex probability rep4": [normalized_complex_probability_values[4]],
                "mean normalized complex probability": [np.mean(normalized_complex_probability_values)],
                "sd normalized complex probability": [np.std(normalized_complex_probability_values)],
                "ttest t_stat": [t_stat],
                "ttest p_value": [p_value],
                "ttest correction": [p_value * len(complexes_hpi)]
            })

            hpi_results = pd.concat([hpi_results, complex_hpi_df], ignore_index=True)

        hpi_results.to_excel(path_results + f"mean_likelihoods_hpis/hpi{hpi}.csv", index=False)
    return

for j in range(1,9):
    paths_rank = [f"results/complex_counts_rep{i}.csv" for i in range(total_reps)]
    calculate_complex_likelihood(paths_rank, f"/home/s391697/f2/complex_likelihood/results/rank{j}/", j, False, False)