# BrainEigenmodes_EDRLR

Code related to the preprint titled "[Beyond cortical geometry: brain dynamics shaped by long-range connections](https://www.biorxiv.org/content/10.1101/2024.04.09.588757v1.abstract)"

## Folder descriptions

1. `Connectome_Derivation/`: folder containing data and codes to derive the various anatomical priors
2. `data/`: folder containing various outputs for the different steps of the pipleine - decomposition, projections, dataset analysis and plotting
3. `results/`: folder containing various outputs to reproduce the Figure 2 and 3 of the manuscript

## File descriptions

1. `generate_Figure2_eigenmode_analysis_32k_EDRLR_all_subjects_fmri.m`: MATLAB script to reproduce Figure 2
2. `generate_result1_eigenmode_analysis_32k_EDRLR_all_subjects_tmri.m`: MATLAB script to reproduce parts of Figure 3
3. `generate_Figure3_eigenmode_plotting_32k_connectome_v7_tk14.m`: MATLAB script to reproduce parts of Figure 3

## Installation
Simply download the repository to get started.
In order to generate the manuscript Figures additional files have to be downloaded from the OSF repository: [https://osf.io/asntf/](https://osf.io/3qjp5/).

## Original data
The original empirical data stem from the Human Connectome Project. Refer to the provided link for comprehensive access, licensing, and usage terms.
Codes and parts of the analysis are built upon Pang et al. 2023 publication https://www.nature.com/articles/s41586-023-06098-1

## Compatibility
Codes are tested on MATLAB versions R2023b.

## Citation

[Preprint] J. Vohryzek, M.L. Kringelbach ,G. Deco, Beyond cortical geometry: brain dynamics shaped by long-range connections, bioRxiv (2004) (DOI: [doi.org/10.1101/2024.04.09.588757](https://www.biorxiv.org/content/10.1101/2024.04.09.588757v1.abstract))

## Further details
contact jakub.vohryzek@upf.edu
