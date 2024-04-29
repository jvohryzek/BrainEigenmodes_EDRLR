%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% calucalation of EDRLR eigenmodes for laplacian project
%
% Jakub Vohryzek June 2023
%
%% loading data
% change directory here
repo_dir = '/Users/jakub/Matlab/Collaboration_Deco/project_laplacian/BrainEigenmodes_EDRLR-main';
addpath(genpath(repo_dir))
%% Load surface files for fsaverage5

surface_interest = 'fsLR_32k';
hemisphere = 'lh';
mesh_interest = 'midthickness';
% Load cortex mask
cortex = dlmread(sprintf('Data/template_surfaces_volumes/%s_cortex-%s_mask.txt', surface_interest, hemisphere));
cortex_ind = find(cortex);
% Calculate number of vertices
num_vertices = length(cortex);
%% load EDRLR
load('Connectome_derivation/Connectomes/EDR_LR_derivation_v2.mat', 'EDR_LRE')

%%
num_modes = 200;

size(EDR_LRE)
[eig_vec_temp, eig_val] = calc_network_eigenmode(EDR_LRE, num_modes);
% Bring back medial wall vertices with zero values
eig_vec = zeros(num_vertices, num_modes);
eig_vec(cortex_ind,:) = eig_vec_temp(:,1:num_modes);

save(sprintf('Connectome_derivation/Connectomes/synthetic_EDRLR_version2_eigenmodes-%s_%i.mat', hemisphere, num_modes), 'eig_val', 'eig_vec', '-v7.3')
