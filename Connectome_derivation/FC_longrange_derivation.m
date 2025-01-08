%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calucalation of parcellateed high density long-range exception EDR connnectomee
%%%
%%% Original: Jakub Vohryzek, Universitat Pompeu Fabra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% loading data
% change directory here
repo_dir = '/project_laplacian/BrainEigenmodes_EDRLR-main';
addpath(genpath(repo_dir))

%% Load surface files for fsaverage5

surface_interest = 'fsLR_32k';
hemisphere = 'lh';
mesh_interest = 'midthickness';
EDR_fitting = 0; % 1 - yes, 0 - no
num_sbj = 255;
[vertices, faces] = read_vtk(sprintf('Data/template_surfaces_volumes/%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';

% Load cortex mask
cortex = dlmread(sprintf('Data/template_surfaces/%s_cortex-%s_mask.txt', surface_interest, hemisphere));
cortex_ind = find(cortex);

% Calculate number of vertices
num_vertices = length(cortex);

% Calculate synthetic connectome eigenmodes

hemisphere = 'lh';
surface_to_analyze = surface_midthickness;

triu_ind_vertex= calc_triu_ind(zeros(size(cortex_ind,1), size(cortex_ind,1)));
% calculated parcellated data
% At parcellated level
parc_name = 'Glasser360';
parc = dlmread(sprintf('Data/parcellations/fsLR_32k_%s-%s.txt', parc_name, hemisphere));
num_parcels = length(unique(parc(parc>0)));

% Extract upper triangle indices
triu_ind = calc_triu_ind(zeros(num_parcels, num_parcels));
%% subject list

% Read the subject IDs from the text file
subjectList = importdata('Data/empirical/subject_list_HCP.txt');

%% loading the empirical connectome and parcellating
load('empirical/S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat')
connectome = avgSC_L;
 
% parcellated
connectome_parc = calc_parcellate_matrix(parc(cortex_ind), connectome);
connectome_parc_matrix = connectome_parc./max(max(connectome_parc));
clear avgSC_L connectome
connectome_parc_vec = connectome_parc_matrix(triu_ind);
%connectome_parc_vec = connectome_parc_vec./max(connectome_parc_vec);

%% calculating the connections euclidean distance

surface_dist = squareform(pdist(surface_to_analyze.vertices(cortex_ind,:)));

% parcellated
surface_dist_parc = calc_parcellate_matrix(parc(cortex_ind), surface_dist);
clear surface_dist
surface_dist_parc_vec = surface_dist_parc(triu_ind);

FC_LRE_euc_dist_vec    = surface_dist_parc_vec;
FC_LRE_euc_dist_matrix = surface_dist_parc;


%%
for sbj = 1:num_sbj
    disp(['subject #',num2str(sbj)])

    % load all-subjects= rfMRI timeseries data
    data = load(['/Users/jakub/Datasets/HCP/HCP255/subject_',num2str(subjectList(sbj)),sprintf('_rfMRI_timeseries-%s.mat', hemisphere)]);
    data_to_reconstruct = data.timeseries;
    T = size(data_to_reconstruct, 2);

    % Calculate empirical FC
    data_parc_emp = calc_parcellate(parc, data_to_reconstruct);
    data_parc_emp = calc_normalize_timeseries(data_parc_emp');
    data_parc_emp(isnan(data_parc_emp)) = 0;
    
    FC_emp_temp = data_parc_emp'*data_parc_emp;
    FC_emp_temp = FC_emp_temp/T; % empirical FC
    FCvec_emp = FC_emp_temp(triu_ind); % this is what is being fitted
    FC_emp_temp_HCP_255_subject(sbj,:,:) = FC_emp_temp;
    FCvec_emp_temp_HCP_255_subject(sbj,:) = FCvec_emp;

end
FC_emp_HCP_255_subject_mean = squeeze(mean(FC_emp_temp_HCP_255_subject,1));
FCvec_emp_subject_mean = squeeze(mean(FCvec_emp_temp_HCP_255_subject,1));

%% results
FC_LRE_euc_dist_vec_thresholded = cell(1,4);
FC_LRE_euc_dist_matrix_thresholded = cell(1,4);

th = 40;%
for i=1:4
    FC_LRE_euc_dist_vec_thresholded{i}    = (FC_LRE_euc_dist_vec).*((FC_LRE_euc_dist_vec)>th(i));
    FC_LRE_euc_dist_matrix_thresholded{i} = (FC_LRE_euc_dist_matrix).*((FC_LRE_euc_dist_matrix)>th(i));%% plotting
end
%% Summary of euclidean distance and FCemp

figure
subplot(2,2,1)
imagesc(FC_LRE_euc_dist_matrix);title('Euclidean Distance')
axis square
subplot(2,2,2)
imagesc(FC_LRE_FC_corr_matrix);title('FC emp')
axis square
subplot(2,2,3)
histogram(nonzeros(FC_LRE_euc_dist_vec));
xlabel('euclidean_distance','Interpreter','None')
ylabel('number of connections');axis square
subplot(2,2,4)
histogram(nonzeros(FC_LRE_FC_corr_vec))
xlabel('FC emp','Interpreter','None')
ylabel('number of connections');axis square

%% figures euclidean distance across thresholds

figure
for i =1:4
    subplot(2,2,i)
    imagesc(FC_LRE_euc_dist_matrix_thresholded{i});axis square
    xlabel('euclidean_distance','Interpreter','None')
    ylabel('number of connections');axis square

end

%% figures across thresholds FC emp

figure,
subplot(1,5,1)
imagesc(FC_LRE_FC_corr_matrix);axis square;colorbar

for i=1:4
    th1 = [0, 0.2, 0.5, 0.8]
    subplot(1,5,i+1);
    imagesc(FC_LRE_FC_corr_matrix>th1(i));axis square;colorbar
end