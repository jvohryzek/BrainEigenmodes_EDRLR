# BrainEigenmodes_EDRLR

Code related to the preprint titled "[Human brain dynamics are shaped by rare long-range connections over and above cortical geometry](https://www.pnas.org/doi/10.1073/pnas.2415102122)"

Short Summary

Explaining how structure of the brain gives rise to its emerging dynamics is a primary pursuit in neuroscience. We describe a fundamental anatomical constraint that emphasizes the key role of rare long-range (LR) connections in explaining functional organization of the brain in terms of spontaneous and task-evoked activity. Specifically, this constraint unifies brain geometry and local connectivity through the exponential distance rule while considering the LR exceptions to this local connectivity as derived from the structural connectome. In addition, when using this structural information, we show that the task-evoked brain activity is described by a low-dimensional manifold of several modes, suggesting that less is more for the efficient information processing in the brain.

## Folder descriptions

1. `Connectome_Derivation/`: folder containing data and codes to derive the various anatomical priors
2. `data/`: folder containing various outputs for the different steps of the pipleine - decomposition, projections, dataset analysis and plotting
3. `results/`: folder containing various outputs to reproduce the Figure 2 and 3 of the manuscript

## File descriptions

1. `PNAS_Figure_2_main.m`: MATLAB script to reproduce Figure 2
2. `PNAS_Figure_3_main_part_ABC.m`: MATLAB script to reproduce parts of Figure 3
3. `PNAS_Figure_3_main_part_CD.m`: MATLAB script to reproduce parts of Figure 3

## Installation
Simply download the repository to get started.
In order to generate the manuscript Figures additional files have to be downloaded from the OSF repository: [https://osf.io/asntf/](https://osf.io/3qjp5/).

## Original data
The original empirical data stem from the [Human Connectome Project](https://www.humanconnectome.org/). Refer to the provided link for comprehensive access, licensing, and usage terms.

Codes and parts of the analysis are built upon Pang et al. 2023 publication https://www.nature.com/articles/s41586-023-06098-1

## Compatibility
Codes are tested on MATLAB versions R2023b.

## Citation

[PNAS] J. Vohryzek, Sanz-Perl Y, M.L. Kringelbach ,G. Deco, Human brain dynamics are shaped by rare long-range connections over and above cortical geometry, PNAS (2025) (DOI: [https://doi.org/10.1073/pnas.2415102122](https://www.pnas.org/doi/10.1073/pnas.2415102122))

## Further details
contact: jakub.vohryzek@upf.edu
