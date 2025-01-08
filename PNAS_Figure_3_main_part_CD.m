%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 3 ABC
%%% Jakub Vohryzek, Universitat Pompeu Fabra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loading data
repo_dir = '/BrainEigenmodes_EDRLR_PNAS-main';
addpath(genpath(repo_dir))

%% Load surface files for visualization
hemisphere = 'lh';
num_modes = 200;
num_sbj = 255;
num_tk = 47;
surface_interest = 'fsLR_32k';
mesh_interest = 'midthickness';

[vertices, faces] = read_vtk(sprintf('Data/template_surfaces/%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';
%% Load cortex mask
cortex = dlmread(sprintf('Data/template_surfaces/%s_cortex-%s_mask.txt', surface_interest, hemisphere));
cortex_ind = find(cortex);

%% Calculate number of vertices
num_vertices = length(cortex);
%
data = load('Data/empirical/S255_tfMRI_ALLTASKS_raw_lh.mat');
fieldNames = fieldnames(data.zstat);
%% subject list

% Read the subject IDs from the text file
subjectList = importdata('Data/empirical/subject_list_HCP.txt');

%% Representing the activation maps for tMRI

data = load('Data/empirical/S255_tfMRI_ALLTASKS_raw_lh.mat');
fieldNames = fieldnames(data.zstat);
activation_map_motor_rf_avg         = nanmean(data.zstat.motor_rf_avg,2);
activation_map_relational_match_rel = nanmean(data.zstat.relational_match_rel,2);

%%
representTask = [3, 11, 19, 30, 41, 44, 47];
ith = [5,7];


%% load here
load('Results/long_3T_task/reconstructions_200modes_motor_rf_avg.mat')

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
