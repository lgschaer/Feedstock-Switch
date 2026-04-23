# Feedstock-Switch

## Overview

This repository contains the **data and R scripts** used to reproduce the analyses and figures for the manuscript:

**“Acetate to caproate: metagenomic insights into functional shifts in a methane-arrested anaerobic bioreactor”**  
_Currently in revision at FEMS Microbes._

This is a reproducibility-focused repository providing:

- Raw and processed annotation data  
- Figure-generating R scripts  
- Exported figures used in the manuscript  

---

## Repository Structure

### 📊 Figures

The repository includes the final and intermediate versions of figures used in the manuscript:

- `Feed switch graphical abstract.png`
- `this_is_for_biorender_graph_abs.png`
- `Figure 1A v1.png`, `Figure 1A v2.png`
- `Figure 1AB v1.png`, `Figure 1AB v1.tiff`
- `Figure 1B v3.png`
- `Figure 2.png`
- `Figure 3 v2.png`
- `Figure 4 v2.png`

These correspond to the main figures presented in the paper.

---

### 📈 R Scripts

Each script generates a specific figure or analysis:

- `Figure1_VFAs_RelAb.R`  
  → Relative abundance of volatile fatty acids (Figure 1)

- `Figure2_FnRedundancy.R`  
  → Functional redundancy analysis (Figure 2)

- `Figure3_Phenotype.R`  
  → Phenotypic profiling (Figure 3)

- `Figure4_Bacteriocins.R`  
  → Bacteriocin-related analysis (Figure 4)

- `parsing_EggNOG-mapper_out amino acids SI only.R`
   → Amino a cid pathway analysis (SI only)

These scripts rely on the annotation datasets provided below.

---

### 🧬 Data Files

Files supporting Figure 1: abundance, taxonomy, MAG quality, and volatile fatty acid measurements

- `clean_mag_ids.csv`
- `coverm_relative_abundance_HQ_MQ.csv`
- `taxQC_info.csv`
- `VFA_data_RYAN.csv`

Files supporting Figures 2, 3, and 4: MAG-level information and annotations

- `clean_mag_ids.csv`
- `W_cat.annotations.csv`
- `Y5_Co_cat.annotations.csv`
- `H_cat.annotations.csv.zip`
- `Y5_8_cat.annotations.csv.zip`
- `Y5_9_cat.annotations.csv.zip`
- `metaErg_clean_genes_out_09192025.csv.zip`
- `merged_and_cleaned_eggnog_output_11192025.csv.zip`

These files contain gene annotations and categorical assignments used in downstream analyses.

> Note: Some files are compressed (`.zip`) due to size.

---

## Reproducibility

To reproduce the figures:

1. Download or clone the repository  
2. Unzip all `.zip` data files  
3. Open the relevant `.R` script  
4. Ensure required R packages are installed  
5. Run the script to regenerate the corresponding figure  

Each script is designed to work directly with the provided data files.

---

## Citation

If you use these data or scripts, please cite the associated manuscript:

> Schaerer et al. 2026
> *Acetate to caproate: metagenomic insights into functional shifts in a methane-arrested anaerobic bioreactor*  
> (currently in revision at FEMS Microbes)

---

## Contact

For questions about the data or analysis, please contact the corresponding author via GitHub or laura.schaerer@colostate.edu
