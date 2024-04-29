%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EDR_connectomes_derivation.m
%%%
%%% Original: James Pang, Monash University
%%% updated: Jakub Vohryzek, Universitat Pompeu Fabra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loading data
% change directory here
repo_dir = '/Users/jakub/Matlab/Collaboration_Deco/project_laplacian/BrainEigenmodes_EDRLR-main';
addpath(genpath(repo_dir))
%% Load surface files for fsaverage5

surface_interest = 'fsLR_32k'; % fsLR_32k 
hemisphere = 'lh';
mesh_interest = 'midthickness';

[vertices, faces] = read_vtk(sprintf('Data/template_surfaces_volumes/%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';
% Load cortex mask
cortex = dlmread(sprintf('Data/template_surfaces_volumes/%s_cortex-%s_mask.txt', surface_interest, hemisphere));
cortex_ind = find(cortex);

% Calculate number of vertices
num_vertices = length(cortex);
%% Calculate synthetic EDR connectome eigenmodes

surface_to_analyze = surface_midthickness;

% =========================================================================
% Calculate Euclidean distance of surface vertices without the medial wall                    
% =========================================================================

surface_dist = squareform(pdist(surface_to_analyze.vertices(cortex_ind,:)));

% =========================================================================
%                    Generate synthetic EDR connectome                     
% =========================================================================

% Probability function
Pspace_func = @(scale, distance) exp(-scale*distance);

% Generate pseudorandom numbers to compare with probability function
% random number is set for now for reproducibility
rng(1)
rand_prob = rand(size(surface_dist));
rand_prob = triu(rand_prob,1) + triu(rand_prob,1)';
rand_prob(1:1+size(rand_prob,1):end) = 1;

% Calculate probability
% Scale = 0.120 matches empirical structural connectivity data. But you can 
% change its value accordingly.
scale = 0.120;
Pspace = Pspace_func(scale, surface_dist);
Pspace(1:1+size(Pspace,1):end) = 0;
Pspace = Pspace/max(Pspace(:));

for belkin_imp = 1:2
    
    if belkin_imp == 1 % EDR binary
        % Generate EDR connectome
        connectome = double(rand_prob < Pspace);
        connectome(1:1+size(connectome,1):end) = 0;
    elseif belkin_imp == 2  % EDR continuous
        connectome = Pspace;
        connectome(1:1+size(connectome,1):end) = 0;
    end

    %                           Calculate the modes                            
    % =========================================================================
    
    num_modes = 200;
    [eig_vec_temp, eig_val] = calc_network_eigenmode(connectome, num_modes);
    
    % Bring back medial wall vertices with zero values
    eig_vec = zeros(num_vertices, num_modes);
    eig_vec(cortex_ind,:) = eig_vec_temp(:,1:num_modes);
    if belkin_imp == 1
        save(sprintf('Connectome_derivation/Connectomes/synthetic_EDRconnectome_binary_eigenmodes-%s_%i.mat', hemisphere, num_modes), 'eig_val', 'eig_vec', '-v7.3')
    elseif belkin_imp == 2
        save(sprintf('Connectome_derivation/Connectomes/synthetic_EDRconnectome_weighted_eigenmodes-%s_%i.mat', hemisphere, num_modes), 'eig_val', 'eig_vec', '-v7.3')
    end
end
% =========================================================================
%                      Some visualizations of results                      
% =========================================================================
%%
% 1st to 5th modes with medial wall view
mode_interest = [51:55]; % [1:5] or [51:55] or [101:105] or  [151:155]
surface_to_plot = surface_midthickness;
data_to_plot = eig_vec(:, mode_interest);
medial_wall = find(cortex==0);
with_medial = 1;

fig = draw_surface_bluewhitered_gallery_dull(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial);
fig.Name = 'Multiple EDR connectome eigenmodes with medial wall view';
