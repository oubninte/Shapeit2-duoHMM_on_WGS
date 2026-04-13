
# Shapeit2-duoHMM on WGS: Recombination Inference and Analysis
    
> Comprehensive R and Shell scripts for analyzing the reliability and accuracy of recombination inference using Shapeit2 duoHMM on whole genome sequence (WGS) data.

## 📋 Overview

This repository contains all computational code and analyses from *"The reliability and accuracy of recombination inferred by Shapeit2 duoHMM on whole genome sequence"* study. The project focuses on inferring recombination events from WGS data using the Shapeit2 duoHMM algorithm and comparing results with linkage analysis methods (Merlin).

**Key Features:**
    - Pedigree structure manipulation and family subdivision
    - Recombination inference from Shapeit2 and OmniExpress data
    - Comparative analysis with Merlin linkage analysis results
    - Statistical validation and overlap analysis
    - Comprehensive pedigree quality control and visualization



## 🔧 Main Components

## 1. **Pedigree Management**

### Cutting_LargeFamillies.R : Subdivides large families into manageable sub-families for linkage analysis.

### peds_statistics.R : Computes comprehensive descriptive statistics for pedigree structures.

### 2. **Recombination Analysis**

    ### Recomb_Number_Analyzes_allchr.R
    Analyzes the count distribution of recombination events across all chromosomes.
    
    ### Recomb_Overlap_Analyzes_allchr.R
    Compares recombination events detected between different methods.
    
    ### OverlapConfidence.R
    Calculates confidence metrics for overlapping recombination regions.
    
    ### Recomb_EA-AA_Analyzes2.R
    Stratified analysis by population ancestry (European Ancestry and African Ancestry populations).

### 3. Merlin recombination events

    ### Recomb_Merlin_UncutFam_allchr.R
    Recombination events inférence by Merlin on Omniexpress chip
    
    ### Merlin_Recomb_allchr.sh
    Shell script orchestrating Merlin recombination events inference pipeline.


📧 Contact & Support
For questions regarding this code, please open an issue in the GitHub repository or contact the repository maintainer.

🙏 Acknowledgments
This work utilizes the following tools and packages:

    R Core Team for R statistical computing
    Goncalo Abecasis for Merlin linkage analysis software
    Jonathan Marchini for Shapeit2 haplotype inference
    Package developers (dplyr, xtable, kinship2, FamAgg)

