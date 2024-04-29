# BrainEigenmodes_EDRLR

Code related to the preprint titled "[Beyond cortical geometry: brain dynamics shaped by long-range connections](https://www.biorxiv.org/content/10.1101/2024.04.09.588757v1.abstract)"

The Harmonic Decomposition of Spacetime (HADES) framework that characterises how different harmonic modes defined in space are expressed over time. Here applied to the DMT dataset

## Folder descriptions

1. `data/`: folder containing data for the different steps of the HADES pipleine - decomposition, projections, dataset analysis and plotting
2. `utils/`: folder containing matlab functions for the different steps of the HADES pipleine - decomposition, projections, dataset analysis and plotting
3. `results/`: folder containing various outputs for the different steps of the HADES pipleine - decomposition, projections, dataset analysis and plotting
4. `figures/`: folder containing figures for the different steps of the HADES pipleine - decomposition, projections, dataset analysis and plotting

## File descriptions


1. `p1_HADES_DMT_FMRI_main_projectFH.m`: MATLAB script to project functional harmonics onto the timeseries
2. `p2_HADES_DMT_spatiotemporal_analysis.m`: MATLAB script to calculate the spatio-temporal analysis
3. `p3_HADES_DMT_dynamic_analysis_publication.m`: MATLAB script to calculate dynamic analysis
4. `p4_HADES_DMT_latent_space_analysis_publication.m`: MATLAB script to calculate latent space analysis analysis

## Supplementary file descriptions
1. `s1_HCP_denseFC_2_vertices.m`: MATLAB script to load the dense FC
2. `s2_HADES_basis_denseFC_vertex_on_HCP.m`: MATLAB scripts to run the laplace decomposition on the dense FC
3. `s3_HADES_plotting_basis.m`: MATLAB script to plot the functional harmonics on the cortical surface

## Installation
Simply download the repository to get started.
In order to run the code two additional files (Functional Harmonics (FHs) and projections of FHs to the DMT dataset) have to be downloaded from the OSF repository: https://osf.io/asntf/ you can find further instruction in the results folder.

HADES pipeline [p1-p4].

To run the code p1 uses the FHs projections on the fMRI resting-state data. The OSF repository provided the FHs projections onto the DMT dataset.

For personal use the appropriate dataset FHs projections should be provided in the cifti 64k vertices format.

Inside each code file, you'll find comments and documentation to guide you through usage.
The repository serves as standalone for the HADES method. Please Consult the documentiaton for further guidance

## Downloading data
Due to privacy the DMT data is only provided in terms of the FHs projections. The original data is available upon request from the authors of the experiment.
To derive the Functional Harmonics, the dense functional connectome of the HCP dataset was used and can be accessed here.
Important: Certain parts of generate_paper_figures.m and generate_paper_suppfigures.m rely on this OSF-hosted data. Ensure it's saved in the correct folders for smooth script execution.

## Original data
The original empirical data stem from the Human Connectome Project. Refer to the provided link for comprehensive access, licensing, and usage terms.


## Compatibility
Codes are tested on MATLAB versions R2023b.

## Citation

[Preprint] J. Vohryzek, M.L. Kringelbach ,G. Deco, Beyond cortical geometry: brain dynamics shaped by long-range connections, bioRxiv (2004) (DOI: [doi.org/10.1101/2024.04.09.588757](https://www.biorxiv.org/content/10.1101/2024.04.09.588757v1.abstract))

## Further details
contact jakub.vohryzek@upf.edu
