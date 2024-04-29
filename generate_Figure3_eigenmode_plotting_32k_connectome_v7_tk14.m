%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loading data
% change directory here
repo_dir = '/Users/jakub/Matlab/Collaboration_Deco/project_laplacian/BrainEigenmodes_EDRLR-main';
addpath(genpath(repo_dir))

%% Load surface files for visualization
hemisphere = 'lh';
num_modes = 200;
num_sbj = 255;
num_tk = 47;
surface_interest = 'fsLR_32k';
mesh_interest = 'midthickness';

[vertices, faces] = read_vtk(sprintf('Data/template_surfaces_volumes/%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';
%% Load cortex mask
cortex = dlmread(sprintf('Data/template_surfaces_volumes/%s_cortex-%s_mask.txt', surface_interest, hemisphere));
cortex_ind = find(cortex);

%% Calculate number of vertices
num_vertices = length(cortex);
%
data = load('Data/empirical/S255_tfMRI_ALLTASKS_raw_lh.mat');
fieldNames = fieldnames(data.zstat);
%% subject list

% Read the subject IDs from the text file
subjectList = importdata('Data/empirical/subject_list_HCP.txt');
%%
recon_temp = zeros(size(cortex_ind,1),2,num_sbj,num_modes);
recon_corr_vertex = zeros(2,num_sbj,num_modes);
recon_beta = zeros(2,num_sbj,num_modes, num_modes);
%% Representing the activation maps for tMRI

data = load('Data/empirical/S255_tfMRI_ALLTASKS_raw_lh.mat');
fieldNames = fieldnames(data.zstat);
activation_map_motor_rf_avg         = nanmean(data.zstat.motor_rf_avg,2);
activation_map_relational_match_rel = nanmean(data.zstat.relational_match_rel,2);

%%
representTask = [3, 11, 19, 30, 41, 44, 47];
ith = [5,7];


%% run above or load here
load('Results/reconstructions_200modes_motor_rf_avg.mat')


%% rendering activation_map_relational_match_rel (Figure 3C)
surface_to_plot = surface_midthickness;
data_to_plot = activation_map_relational_match_rel;
medial_wall = find(cortex==0);
with_medial = 1;

fig = draw_surface_bluewhitered_gallery_dull(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial);
fig.Name = 'activation_map_motor_rf_avg';

%% Figure 3C: EDRLR - Reconstrcuting specific brain modes

data_connectome_file = sprintf('Connectome_derivation/Connectomes/synthetic_EDRLR_eigenmodes_fsLR_32k-%s_%i.mat', hemisphere, num_modes);
data_connectome = load(data_connectome_file);
eig_vec_EDRLR = data_connectome.eig_vec;
%
mode_interest =[16]; % [2, 3, 4, 5, 7, 16] % the specific brain modes
surface_to_plot = surface_midthickness;
data_to_plot = eig_vec_EDRLR(:, mode_interest);
medial_wall = find(cortex==0);
with_medial = 1;

fig = draw_surface_bluewhitered_gallery_dull(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial);
fig.Name = 'Multiple EDRLR connectome eigenmodes with medial wall view';
%% Figure 3D:  activation_map_motor_rf_avg 
surface_to_plot = surface_midthickness;
data_to_plot = activation_map_motor_rf_avg;
medial_wall = find(cortex==0);
with_medial = 1;

fig = draw_surface_bluewhitered_gallery_dull(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial);
fig.Name = 'activation_map_motor_rf_avg';
%% Figure 3D: Geometry - Reconstrcuting activation_map_motor_rf_avg 
mode_interest = [20,15,10,5];
surface_to_plot = surface_midthickness;
data_to_plot = recon_temp_cortex_v5(:,mode_interest);
medial_wall = find(cortex==0);
with_medial = 1;

fig = draw_surface_bluewhitered_gallery_dull(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial);
fig.Name = 'Geometry - activation_map_motor_rf_avg';
%% Figure 3D: EDRLR - Reconstrcuting activation_map_motor_rf_avg 
mode_interest = [20, 15, 10, 5];
surface_to_plot = surface_midthickness;
data_to_plot = recon_temp_cortex_v7(:,mode_interest);
medial_wall = find(cortex==0);
with_medial = 1;

fig = draw_surface_bluewhitered_gallery_dull(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial);
fig.Name = 'EDRLR - activation_map_motor_rf_avg';
%%
recon_beta_v5  = squeeze(mean(recon_beta(1,:,:,:),2));
% Define the number of categories
numCategories = 20;

% Create data for the spider plot 
data = diag(abs(recon_beta_v5));
data = data(1:20)';
% Create angles for each category
angles = linspace(0, 2*pi, numCategories + 1);
angles(end) = [];  % Remove the last angle to close the plot

% Create a figure
figure;

% Plot the spider plot
polarplot([angles, 0], [data, data(1)], 'k-','LineWidth',2);  % Repeat the first data point to create a closed loop

% Customize the plot
title('Geometry - Motor rf avg');
thetaticks(angles.*(180/pi))
thetaticklabels({ 'Mode 0','Mode 1', 'Mode 2', 'Mode 3', 'Mode 4', 'Mode 5', ...
    'Mode 6', 'Mode 7', 'Mode 8', 'Mode 9', 'Mode 10', ...
    'Mode 11', 'Mode 12', 'Mode 13', 'Mode 14', 'Mode 15', ...
    'Mode 16', 'Mode 17', 'Mode 18', 'Mode 19'});
% Define the number of categories
numCategories = 20;
%%
recon_beta_v7  = squeeze(mean(recon_beta(2,:,:,:),2));

% Create data for the spider plot 
data = diag(abs(recon_beta_v7));
data = data(1:20)';
% Create angles for each category
angles = linspace(0, 2*pi, numCategories + 1);
angles(end) = [];  % Remove the last angle to close the plot

% Create a figure
figure;

% Plot the spider plot
polarplot([angles, 0], [data, data(1)], 'k-','LineWidth',2);  % Repeat the first data point to create a closed loop

% Customize the plot
title('EDRLR - Motor rf avg');
thetaticks(angles.*(180/pi))
thetaticklabels({ 'Mode 0','Mode 1', 'Mode 2', 'Mode 3', 'Mode 4', 'Mode 5', ...
    'Mode 6', 'Mode 7', 'Mode 8', 'Mode 9', 'Mode 10', ...
    'Mode 11', 'Mode 12', 'Mode 13', 'Mode 14', 'Mode 15', ...
    'Mode 16', 'Mode 17', 'Mode 18', 'Mode 19'});

%% %%%%%%%%%%%%%%%%%%%%%%% rendering modes %%%%%%%%%%%%%%%%%%%%%%%
%% rendering modes EDRLR

surface_to_plot = surface_midthickness;
data_to_plot = activation_map_motor_rf_avg;
medial_wall = find(cortex==0);
with_medial = 1;
mode_interest = [2, 3, 4, 5, 7, 16];

fig = draw_surface_bluewhitered_gallery_dull(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial);
fig.Name = 'activation_map_motor_rf_avg';
surface_to_plot = surface_midthickness;
data_to_plot = eig_vec_EDRLR(:, mode_interest);
medial_wall = find(cortex==0);
with_medial = 1;

fig = draw_surface_bluewhitered_gallery_dull(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial);
fig.Name = 'Multiple EDRLR connectome eigenmodes with medial wall view';
%% rendering modes Geometry
data_connectome_file = sprintf('/Users/jakub/Matlab/Collaboration_Deco/project_laplacian/output/synthetic_pang_geometry_eigenmodes_fsLR_32k-%s_%i.mat', hemisphere, num_modes);
data_connectome = load(data_connectome_file);
eig_vec_Geometry = data_connectome.eig_vec;

surface_to_plot = surface_midthickness;
data_to_plot = activation_map_motor_rf_avg;
medial_wall = find(cortex==0);
with_medial = 1;
mode_interest = [2, 3, 4, 5, 7, 16];

fig = draw_surface_bluewhitered_gallery_dull(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial);
fig.Name = 'activation_map_motor_rf_avg';
surface_to_plot = surface_midthickness;
data_to_plot = eig_vec_Geometry(:, mode_interest);
medial_wall = find(cortex==0);
with_medial = 1;

fig = draw_surface_bluewhitered_gallery_dull(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial);
fig.Name = 'Multiple Geometry connectome eigenmodes with medial wall view';
