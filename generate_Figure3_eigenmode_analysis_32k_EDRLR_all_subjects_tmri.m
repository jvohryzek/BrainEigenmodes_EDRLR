%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% generate_result1_eigenmode_analysis_32k_EDRLR_all_subjects_tmri.m
%%%
%%% MATLAB script to demonstrate how to use surface eigenmodes to analyze 
%%% fMRI data. In particular, the script demonstrates how to
%%% (1) reconstruct a task fMRI spatial map,
%%% (2) reconstruct a resting-state fMRI spatiotemporal map and functional
%%%     connectivity (FC) matrix, and
%%% (3) calculate the eigenmode-based power spectral content of a spatial map
%%%
%%% NOTE 1: The script can also be used to analyze fMRI data using other
%%%         types of surface eigenmodes (e.g., connectome eigenmodes). 
%%%         Just change the eigenmodes variable below. However, make sure
%%%         that the variable is an array of size
%%%         [number of vertices x number of modes]. 
%%% NOTE 2: Current demo uses 50 modes. For a proper analysis, we advise 
%%%         using between 100 to 200 modes. 200 template geometric 
%%%         eigenmodes are provided in data/template_eigenmodes.
%%%
%%% Original: James Pang, Monash University (demo_eigenmode_analysis.m)
%%% updated: Jakub Vohryzek, Universitat Pompeu Fabra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% loading data
% change directory here
repo_dir = '/Users/jakub/Matlab/Collaboration_Deco/project_laplacian/BrainEigenmodes_EDRLR-main';
addpath(genpath(repo_dir))

%% Load surface files for visualization
hemisphere = 'lh';
num_modes = 200;
num_sbj = 255;
num_tk = 47;
MRI_type = 'tMRI' % 'fMRI' or 'tMRI'
typeRecon = 'original' % 'parcellated' or 'original'
surface_interest = 'fsLR_32k';
mesh_interest = 'midthickness';

[vertices, faces] = read_vtk(sprintf('Data/template_surfaces_volumes/%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';
% Load cortex mask
cortex = dlmread(sprintf('Data/template_surfaces_volumes/%s_cortex-%s_mask.txt', surface_interest, hemisphere));
cortex_ind = find(cortex);

% Calculate number of vertices
num_vertices = length(cortex);
%% parcellated template
% At parcellated level
parc_name = 'Glasser360';
parc = dlmread(sprintf('Data/parcellations/fsLR_32k_%s-%s.txt', parc_name, hemisphere));
num_parcels = length(unique(parc(parc>0)));
%%
data = load('Data/empirical/S255_tfMRI_ALLTASKS_raw_lh.mat');
fieldNames = fieldnames(data.zstat);
%% subject list

% Read the subject IDs from the text file
subjectList = importdata('Data/empirical/subject_list_HCP.txt');

%% 
% =========================================================================
%                    Load eigenmodes and empirical data                    
% =========================================================================

parfor_run = 0;
if parfor_run == 1
for i = [1,2,3,4] 
    data_connectome_file = '';

    %% extract connectomes
    connectome_version = i;
    if connectome_version == 1
        data_connectome_file = sprintf('Connectome_derivation/Connectomes/synthetic_EDRconnectome_binary_eigenmodes-%s_%s_%i.mat',surface_interest, hemisphere, num_modes);
    elseif connectome_version == 2
        data_connectome_file = sprintf('Connectome_derivation/Connectomes/synthetic_EDRconnectome_weighted_eigenmodes-%s_%s_%i.mat',surface_interest, hemisphere, num_modes);
    elseif connectome_version == 3
        data_connectome_file = sprintf('Connectome_derivation/Connectomes/synthetic_pang_geometry_eigenmodes_fsLR_32k-%s_%i.mat', 'lh', 200);
    elseif connectome_version == 4
        data_connectome_file = sprintf('Connectome_derivation/Connectomes/synthetic_EDRLR_version1_eigenmodes-%s_%i.mat', hemisphere, num_modes);
    end
    data_connectome = load(data_connectome_file);
    eigenmodes = data_connectome.eig_vec;
    %% pre-allocation
    recon_corr_vertex_tk_sbj_version1 = zeros(num_tk,num_sbj, num_modes);
    recon_corr_vertex_tk_sbj_version2 = zeros(num_tk,num_sbj, num_modes);
    recon_corr_vertex_tk_sbj_version3 = zeros(num_tk,num_sbj, num_modes);
    recon_corr_vertex_tk_sbj_version4 = zeros(num_tk,num_sbj, num_modes);

    recon_corr_parc_tk_sbj_version1 = zeros(num_tk,num_sbj, num_modes);
    recon_corr_parc_tk_sbj_version2 = zeros(num_tk,num_sbj, num_modes);
    recon_corr_parc_tk_sbj_version3 = zeros(num_tk,num_sbj, num_modes);
    recon_corr_parc_tk_sbj_version4 = zeros(num_tk,num_sbj, num_modes);

    % Load example single-subject tfMRI z-stat data

    data = load('/Data/empirical/S255_tfMRI_ALLTASKS_raw_lh.mat');
    fieldNames = fieldnames(data.zstat);
    for tk = 1:numel(fieldNames)
        for sbj = 1:255
            fieldName = fieldNames{tk};
            disp(['connectome v',num2str(i),' task ',fieldName,' subject #',num2str(sbj)])
            data_to_reconstruct = data.zstat.(fieldName)(:,sbj);
    
            % =========================================================================
            % Calculate reconstruction beta coefficients using 1 to num_modes eigenmodes
            % =========================================================================
            
            recon_beta = zeros(num_modes, num_modes);
            for mode = 1:num_modes
                basis = eigenmodes(cortex_ind, 1:mode);        
                recon_beta(1:mode,mode) = calc_eigendecomposition(data_to_reconstruct(cortex_ind), basis, 'matrix');
            end
            %
            % =========================================================================
            %     Calculate reconstruction accuracy using 1 to num_modes eigenmodes    
            % =========================================================================
            
            % reconstruction accuracy = correlation of empirical and reconstructed data
            
            % At vertex level
            recon_corr_vertex = zeros(1, num_modes);               
            for mode = 1:num_modes
                recon_temp = eigenmodes(cortex_ind, 1:mode)*recon_beta(1:mode,mode);
                recon_corr_vertex(mode) = corr(data_to_reconstruct(cortex_ind), recon_temp);
            end
            
            % At parcellated level
            parc_name = 'Glasser360';
            parc = dlmread(sprintf('Data/parcellations/fsLR_32k_%s-%s.txt', parc_name, hemisphere));
            
            recon_corr_parc = zeros(1, num_modes);               
            for mode = 1:num_modes
                recon_temp = eigenmodes(:, 1:mode)*recon_beta(1:mode,mode);  
                recon_corr_parc(mode) = corr(calc_parcellate(parc, data_to_reconstruct), calc_parcellate(parc, recon_temp));
            end

            if connectome_version == 1
                recon_corr_vertex_tk_sbj_version1(tk,sbj,:) = recon_corr_vertex;
                recon_corr_parc_tk_sbj_version1(tk,sbj,:) = recon_corr_parc;
        
            elseif connectome_version == 2
                recon_corr_vertex_tk_sbj_version2(tk,sbj,:) = recon_corr_vertex;
                recon_corr_parc_tk_sbj_version2(tk,sbj,:) = recon_corr_parc;
            
            elseif connectome_version == 3
                recon_corr_vertex_tk_sbj_version3(tk,sbj,:) = recon_corr_vertex;
                recon_corr_parc_tk_sbj_version3(tk,sbj,:) = recon_corr_parc;
            
            elseif connectome_version == 4
                recon_corr_vertex_tk_sbj_version4(tk,sbj,:) = recon_corr_vertex;
                recon_corr_parc_tk_sbj_version4(tk,sbj,:) = recon_corr_parc;
            end
        end
    end

    if connectome_version == 1
        recon_corr_vertex_all_version1 = squeeze(mean(recon_corr_vertex_tk_sbj_version1,2));
        recon_corr_parc_all_version1 = squeeze(mean(recon_corr_parc_tk_sbj_version1,2));
        saveFileName = 'Results/tMRI_parcellated_reconstruction_all_subjects_EDRbinary';
        parsave(saveFileName, recon_corr_vertex_all_version1, recon_corr_parc_all_version1,...
        recon_corr_vertex_tk_sbj_version1,recon_corr_parc_tk_sbj_version1)
    elseif connectome_version == 2
        recon_corr_vertex_all_version2 = squeeze(mean(recon_corr_vertex_tk_sbj_version2,2));
        recon_corr_parc_all_version2 = squeeze(mean(recon_corr_parc_tk_sbj_version2,2));
        saveFileName = 'Results/tMRI_parcellated_reconstruction_all_subjects_EDRcontinuous';
        parsave(saveFileName, recon_corr_vertex_all_version2, recon_corr_parc_all_version2,...
        recon_corr_vertex_tk_sbj_version2,recon_corr_parc_tk_sbj_version2)
    elseif connectome_version == 3
        recon_corr_vertex_all_version3 = squeeze(mean(recon_corr_vertex_tk_sbj_version3,2));
        recon_corr_parc_all_version3 = squeeze(mean(recon_corr_parc_tk_sbj_version3,2));
        saveFileName = 'Results/tMRI_parcellated_reconstruction_all_subjects_Geometry';
        parsave(saveFileName, recon_corr_vertex_all_version3, recon_corr_parc_all_version3,...
        recon_corr_vertex_tk_sbj_version3,recon_corr_parc_tk_sbj_version3)
    elseif connectome_version == 4 
        recon_corr_vertex_all_version4 = squeeze(mean(recon_corr_vertex_tk_sbj_version4,2));
        recon_corr_parc_all_version4 = squeeze(mean(recon_corr_parc_tk_sbj_version4,2));
        saveFileName = 'Results/tMRI_parcellated_reconstruction_all_subjects_EDRLR';
        parsave(saveFileName, recon_corr_vertex_all_version4, recon_corr_parc_all_version4,...
        recon_corr_vertex_tk_sbj_version4,recon_corr_parc_tk_sbj_version4)
    end
end
elseif parfor_run ==0
    load('Results/tMRI_parcellated_reconstruction_all_subjects_EDRbinary.mat')
    recon_corr_parc_version1   = squeeze(nanmean(recon_corr_parc_tk_sbj_version,2))';
    recon_corr_vertex_version1 = squeeze(nanmean(recon_corr_vertex_tk_sbj_version,2))';
    load('Results/tMRI_parcellated_reconstruction_all_subjects_EDRcontinuous.mat')
    recon_corr_parc_version2   = squeeze(nanmean(recon_corr_parc_tk_sbj_version,2))';
    recon_corr_vertex_version2 = squeeze(nanmean(recon_corr_vertex_tk_sbj_version,2))';
    load('Results/tMRI_parcellated_reconstruction_all_subjects_Geometry.mat')
    recon_corr_parc_version3   = squeeze(nanmean(recon_corr_parc_tk_sbj_version,2))';
    recon_corr_vertex_version3 = squeeze(nanmean(recon_corr_vertex_tk_sbj_version,2))';
    load('Results/tMRI_parcellated_reconstruction_all_subjects_EDRLR.mat')
    recon_corr_parc_version4   = squeeze(nanmean(recon_corr_parc_tk_sbj_version,2))';
    recon_corr_vertex_version4 = squeeze(nanmean(recon_corr_vertex_tk_sbj_version,2))';
end

%% FIGURE A
%% representative tasks
representTask = [3, 11, 19, 30, 41, 44, 47]; % the 7 representative tasks
th_corr = 0.6
nMod = 200;
figure
subplot(2,4,2)
plot(1:nMod, recon_corr_parc_version1(1:nMod,representTask), 'b-', 'linewidth', 2);grid on; title('EDR binary');hold on;xlim([1,nMod])
axis square;ylim([0, 1])
subplot(2,4,3);axis square
plot(1:nMod, recon_corr_parc_version2(1:nMod,representTask), 'g-', 'linewidth', 2');grid on; title('EDR weighted');hold on;xlim([1,nMod])
axis square;ylim([0, 1])
subplot(2,4,1);axis square
plot(1:nMod, recon_corr_parc_version3(1:nMod,representTask), 'm-', 'linewidth', 2);grid on;ylim([0,1]); title('Geometry');hold on;xlim([1,nMod])
axis square;ylim([0, 1])
subplot(2,4,4);
plot(1:nMod, recon_corr_parc_version4(1:nMod,representTask), 'k-', 'linewidth', 2);grid on; title('EDR+LR');xlim([1,nMod])
axis square;ylim([0, 1])
subplot(2,4,6);
plot((2:nMod), diff([zeros(1,7); recon_corr_parc_version1(2:nMod,representTask)]), 'b-', 'linewidth', 2);grid on; title('EDR binary');ylim([0, .1]);xlim([2,nMod])
subplot(2,4,7);
plot((2:nMod), diff([zeros(1,7); recon_corr_parc_version2(2:nMod,representTask)]), 'g-', 'linewidth', 2);grid on; title('EDR weighted');ylim([0, .1]);xlim([2,nMod])
subplot(2,4,5);
plot((2:nMod), diff([zeros(1,7); recon_corr_parc_version3(2:nMod,representTask)]), 'm-', 'linewidth', 2);grid on; title('Geometry');ylim([0, .1]);xlim([2,nMod])
subplot(2,4,8);
plot((2:nMod), diff([zeros(1,7); recon_corr_parc_version4(2:nMod,representTask)]), 'k-', 'linewidth', 2);grid on; title('EDR+LR');ylim([0, .1]);xlim([2,nMod])

%% FIGURE B
figure('Name', 'tfMRI reconstruction - accuracy difference');
th_mode = 20;
customYTicks = 1:47; % [1, 3, 5, 7, 10];
customYTickLabels = fieldNames;
% Set the custom Y-axis ticks and labels
set(gca, 'YTick', customYTicks, 'YTickLabel', customYTickLabels,'TickLabelInterpreter','none');
subplot(2,3,1)
imagesc((recon_corr_parc_version3(1:th_mode,:) - recon_corr_parc_version1(1:th_mode,:))',[-0.25 0.25]);title('EDR binary')
subplot(2,3,2)
imagesc((recon_corr_parc_version3(1:th_mode,:) - recon_corr_parc_version2(1:th_mode,:))',[-0.25 0.25]);title('EDR weighted')
subplot(2,3,3)
imagesc((recon_corr_parc_version3(1:th_mode,:) - recon_corr_parc_version4(1:th_mode,:))',[-0.25 0.25]);title('EDR+LR')

colorbar('location','eastoutside')
colormap(redblue)
subplot(2,3,4)
bar(mean((recon_corr_parc_version3(1:th_mode,:) - recon_corr_parc_version1(1:th_mode,:))'),'k','LineWidth',2);title('EDR binary');ylim([-0.1 0.1]);
subplot(2,3,5)
bar(mean((recon_corr_parc_version3(1:th_mode,:) - recon_corr_parc_version2(1:th_mode,:))'),'k','LineWidth',2);title('EDR weighted');ylim([-0.1 0.1]);
subplot(2,3,6)
bar(mean((recon_corr_parc_version3(1:th_mode,:) - recon_corr_parc_version4(1:th_mode,:))'),'k','LineWidth',2);title('EDR+LR');ylim([-0.1 0.1]);
%% FIGURE C
% for 1 tasks differences
figure
i=7;
bar(diff([0; recon_corr_parc_version4(2:th_mode,representTask(i))]), 'k', 'linewidth', 2);grid on; title('EDR+LR');ylim([0, .25]);xlim([0.5,th_mode-0.5])
title(['Mode contribution to reconstruction for task ', customYTickLabels{representTask(i)}],'Interpreter', 'none');ylabel('FC correlation contribution'); xlabel('Modes')
%% info for FIGURE D
customYTickLabels{14}
recon_corr_parc_version3([5,10,15,20],14)
recon_corr_parc_version4([5,10,15,20],14)
%% plotting the reconstruction error for both binary and weighted
% laplacian modes
% figure('Name', 'tfMRI reconstruction - accuracy');
% subplot(1,4,1)
% plot(1:num_modes, recon_corr_vertex_version1, 'b--', 'linewidth', 2); title('EDR binary');hold on;ylim([0 1])
% plot(1:num_modes, recon_corr_parc_version1, 'b-', 'linewidth', 2);grid on
% axis square
% subplot(1,4,2);axis square
% plot(1:num_modes, recon_corr_vertex_version2, 'g--', 'linewidth', 2); title('EDR weighted');hold on;ylim([0 1])
% plot(1:num_modes, recon_corr_parc_version2, 'g-', 'linewidth', 2');grid on
% axis square
% subplot(1,4,3);axis square
% plot(1:num_modes, recon_corr_vertex_version3, 'm--', 'linewidth', 2); title('Pang et al. 2023');hold on
% plot(1:num_modes, recon_corr_parc_version3, 'm-', 'linewidth', 2);grid on;ylim([0 1])
% axis square
% subplot(1,4,4);axis square
% plot(1:num_modes, recon_corr_vertex_version4, 'k--', 'linewidth', 2); title('EDR+LR');hold on;ylim([0 1])
% plot(1:num_modes, recon_corr_parc_version4, 'k-', 'linewidth', 2);grid on
% axis square


