# Master Thesis Project  
**Development of an Image Analysis and Modelling Framework for Influenza Genome Assembly**

Author: Jannis Witte  
University: Julius-Maximilians-Universität Würzburg  
Program: M.Sc. Biowissenschaften  
Supervisors: Dr. Markus Ankenbrand, Prof. Dr. Sabine Fischer  
Submission: 14.07.2025  

---

## Overview
This repository contains the code and analysis pipeline developed during my Master's thesis *"Development of an Image Analysis and Modelling Framework for Influenza Genome Assembly"*.  

The goal was to analyze in-situ sequencing (ISS) data of influenza-infected MDCK cells to understand the mechanisms of viral ribonucleoprotein (vRNP) genome assembly.  
A custom pipeline was implemented to preprocess fluorescence microscopy images, detect RNA spots, decode signals, and model genome assembly pathways.

---

## Dataset
- **Biological system:** MDCK cells infected with influenza H1N1 strain *A/Puerto Rico/8/34 (PR8)*  
- **Imaging:** Epifluorescence microscopy (Leica DFC9000GT camera, XYZ motorized stage, hardware autofocus)  
- **Staining:** Direct padlock probing (PLP) with barcoded probes for all 8 vRNP segments  
- **Acquisition:**  
  - 0–8 hours post infection (hpi), in hourly steps  
  - 9–12 fields of view (FOVs) per timepoint, partially overlapping  
  - Z-stacks projected to 2D via Maximum Intensity Projection (MIP)  
  - Processed with Instant Computational Clearing (ICC)  

---

## Pipeline
The image analysis pipeline was developed in **Python 3.9** (Ubuntu 24.04.2) and builds on the [Starfish 0.3.2](https://github.com/spacetx/starfish) library with custom modifications.  

**Main steps:**
1. **Image Registration** – Align sequencing rounds to correct shifts.  
2. **Background Filtering** – White Top-Hat filtering to reduce noise.  
3. **Normalization** – Percentile-based intensity scaling.  
4. **Masking** – Separation of high- and low-density signal regions.  
5. **Decoding**  
   - Spot-based decoding (SBD) → more accurate, especially in non-dense regions  
   - Pixel-based decoding (PBD) → less consistent, evaluated but not used for final analysis  
6. **Nucleus Segmentation** – DAPI channel segmentation of nuclei.  

The pipeline outputs decoded RNA spots/segments as the basis for downstream analyses.  

---

## Downstream Analyses
Using the decoded spots and segment assignments produced by the pipeline, the following analyses were performed:  

- **Complex Counts** – quantification of single segments and multi-segment complexes.  
- **Complex Probability Calculation** – estimation of expected vRNP combinations from observed segment abundances.  
- **Assembly Modelling** – simulation of genome assembly pathways under a monomeric model.  

---

## Software & Dependencies
The project uses Python with the following main packages:

- `numpy`, `pandas`, `scipy`, `matplotlib`, `scikit-image`, `scikit-learn`  
- `starfish` (custom fork: [starfish_jannis](https://github.com/janniswi00/starfish_jannis))  
- `napari` (for manual annotation with custom plugin: [np-plugin-manual-spot-annotator](https://github.com/janniswi00/np-plugin-manual-spot-annotator))  

Installation via **mamba/conda** is recommended.