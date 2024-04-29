%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% demo_eigenmode_analysis.m
%%% Original: James Pang, Monash University, 
%%% updated: Jakub Vohryzek, Universitat Pompeu Fabra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% loading data
% change directory here
repo_dir = '/Users/jakub/Matlab/Collaboration_Deco/project_laplacian/BrainEigenmodes_EDRLR-main';
addpath(genpath(repo_dir))

%% Load surface files for visualization
hemisphere = 'lh';
num_modes = 200;
num_sbj = 255 % 50% 255;
MRI_type = 'fMRI'; % 'fMRI' or 'tMRI'
surface_interest = 'fsLR_32k';
mesh_interest = 'midthickness';

% load surface
[vertices, faces] = read_vtk(sprintf('Data/template_surfaces/%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';
% Load cortex mask
cortex = dlmread(sprintf('Data/template_surfaces/%s_cortex-%s_mask.txt', surface_interest, hemisphere));
cortex_ind = find(cortex);

% Calculate number of vertices
num_vertices = length(cortex);

%% parcellated template
% At parcellated level
parc_name = 'Glasser360';
parc = dlmread(sprintf('data/parcellations/fsLR_32k_%s-%s.txt', parc_name, hemisphere));
num_parcels = length(unique(parc(parc>0)));

%% subject list

% Read the subject IDs from the text file
subjectList = importdata('Data/empirical/subject_list_HCP.txt');
%% Extract upper triangle indices
triu_ind = calc_triu_ind(zeros(num_parcels, num_parcels));
        
%% functional connectome long range
%% version 1 and 2 of the SC exceptions
load('Connectome_derivation/Euclidean_distance_LRE.mat')

th_dis = [40, 60, 80, 100];

%%
% =========================================================================
% Calculate reconstruction beta coefficients using 1 to num_modes eigenmodes
% =========================================================================

parLoop = 0;
if parLoop == 1
parfor l = 1:4
% for l = 1:4 

    idx = [1,2,3,4];
    i = idx(l);
    data_connectome_file = '';

    %% extract connectomes
    connectome_version = i;
    if connectome_version == 1
        data_connectome_file = sprintf('Connectome_derivation/Connectomes/synthetic_EDRbinary_eigenmodes-%s-%s_%i.mat',surface_interest, hemisphere, num_modes);
    elseif connectome_version == 2
        data_connectome_file = sprintf('Connectome_derivation/Connectomes/synthetic_EDRcontinuous_eigenmodes-%s-%s_%i.mat',surface_interest, hemisphere, num_modes);
    elseif connectome_version == 3
        data_connectome_file = sprintf('Connectome_derivation/Connectomes/synthetic_Geometry_eigenmodes_%s-%s_%i.mat',surface_interest, hemisphere, num_modes);
    elseif connectome_version == 4
        data_connectome_file = sprintf('Connectome_derivation/Connectomes/synthetic_EDRLR_version2_eigenmodes_%s-%s_%i.mat',surface_interest, hemisphere, num_modes);
    end
    data_connectome = load(data_connectome_file);
    eigenmodes = data_connectome.eig_vec;
    % =========================================================================
    %                    Load empirical data                    
    % =========================================================================
    %% pre-allocation

    recon_corr_parc_version1 = zeros(num_sbj, num_modes);
    recon_corr_parc_version2 = zeros(num_sbj, num_modes);
    recon_corr_parc_version3 = zeros(num_sbj, num_modes);
    recon_corr_parc_version4 = zeros(num_sbj, num_modes);

    recon_mse_parc_version1 = zeros(num_sbj, num_modes);
    recon_mse_parc_version2 = zeros(num_sbj, num_modes);
    recon_mse_parc_version3 = zeros(num_sbj, num_modes);
    recon_mse_parc_version4 = zeros(num_sbj, num_modes);

    recon_corr_parc_SC_exceptions_version1 = zeros(num_sbj,size(th_dis,2), num_modes);
    recon_corr_parc_SC_exceptions_version2 = zeros(num_sbj,size(th_dis,2), num_modes);
    recon_corr_parc_SC_exceptions_version3 = zeros(num_sbj,size(th_dis,2), num_modes);
    recon_corr_parc_SC_exceptions_version4 = zeros(num_sbj,size(th_dis,2), num_modes);

    recon_mse_parc_SC_exceptions_version1 = zeros(num_sbj,size(th_dis,2), num_modes);
    recon_mse_parc_SC_exceptions_version2 = zeros(num_sbj,size(th_dis,2), num_modes);
    recon_mse_parc_SC_exceptions_version3 = zeros(num_sbj,size(th_dis,2), num_modes);
    recon_mse_parc_SC_exceptions_version4 = zeros(num_sbj,size(th_dis,2), num_modes);
    %%
    for sbj = 1:num_sbj
        disp(['connectome v',num2str(i),' subject #',num2str(sbj)])
    
        % load all-subjects= rfMIR timeseries data
        data = load(['/Users/jakub/Datasets/HCP/HCP255/subject_',num2str(subjectList(sbj)),sprintf('_rfMRI_timeseries-%s.mat', hemisphere)]);
        data_to_reconstruct = data.timeseries;
        T = size(data_to_reconstruct, 2);
        recon_beta = zeros(num_modes, T, num_modes);
        for mode = 1:num_modes
            basis = eigenmodes(cortex_ind, 1:mode);     
            recon_beta(1:mode,:,mode) = calc_eigendecomposition(data_to_reconstruct(cortex_ind,:), basis, 'matrix');
        end

        % =========================================================================
        %     Calculate reconstruction accuracy using 1 to num_modes eigenmodes    
        % =========================================================================
        
        % reconstruction accuracy = correlation of empirical and reconstructed data

        % Calculate empirical FC
        data_parc_emp = calc_parcellate(parc, data_to_reconstruct);
        data_parc_emp = calc_normalize_timeseries(data_parc_emp');
        data_parc_emp(isnan(data_parc_emp)) = 0;
        
        FC_emp_temp = data_parc_emp'*data_parc_emp;
        FC_emp_temp = FC_emp_temp/T; % empirical FC
        FCvec_emp = FC_emp_temp(triu_ind); % this is what is being fitted
    
        % Calculate reconstructed FC and accuracy (slow to run with more modes)
        FCvec_recon = zeros(length(triu_ind), num_modes);
        recon_corr_parc = zeros(1, num_modes);    
        recon_mse_parc = zeros(1, num_modes);               
        recon_corr_parc_SC_exceptions = zeros(size(th_dis,2),num_modes);
        recon_mse_parc_SC_exceptions = zeros(size(th_dis,2),num_modes);

        for mode = 1:num_modes 
            recon_temp = eigenmodes(:, 1:mode)*squeeze(recon_beta(1:mode,:,mode)); % no need to exclude the midline here because the callc_parcealtew takes care of it. only for vertex resuts
    
            % parcellated
            data_parc_recon = calc_parcellate(parc, recon_temp);
            data_parc_recon = calc_normalize_timeseries(data_parc_recon');
            data_parc_recon(isnan(data_parc_recon)) = 0;
        
            FC_recon_temp = data_parc_recon'*data_parc_recon;
            FC_recon_temp = FC_recon_temp/T;
        
            FCvec_recon(:,mode) = FC_recon_temp(triu_ind);
                   
            recon_corr_parc(mode) = corr(FCvec_emp, FCvec_recon(:,mode));
            recon_mse_parc(mode) = immse(FCvec_emp, FCvec_recon(:,mode));

            % functional long-range exceptions
            for t = 1:4
                %SC_exceptions_matrix = ((SC_exceptions_euc_dist_part1_version2_matrix>th_dis(t)).*(SC_exceptions_FC_corr_part2_version2_matrix>0.5)); 
                FC_LRE_matrix = ((FC_LRE_euc_dist_matrix>th_dis(t)).*(FC_emp_HCP_255_subject_mean>0.5)); 
                a = FC_LRE_matrix.*FC_emp_temp;
                b = FC_LRE_matrix.*FC_recon_temp;
                recon_corr_parc_SC_exceptions(t,mode)  = corr(nonzeros(a(triu_ind)), nonzeros(b(triu_ind)));
                recon_mse_parc_SC_exceptions(t,mode)   = immse(nonzeros(a(triu_ind)), nonzeros(b(triu_ind)));

            end
        end
    
        if connectome_version == 1
            recon_corr_parc_version1(sbj,:) = recon_corr_parc;
            recon_mse_parc_version1(sbj,:)  = recon_mse_parc;        
            recon_corr_parc_SC_exceptions_version1(sbj,:,:) = recon_corr_parc_SC_exceptions;
            recon_mse_parc_SC_exceptions_version1(sbj,:,:) = recon_mse_parc_SC_exceptions;

        elseif connectome_version == 2
            recon_corr_parc_version2(sbj,:) = recon_corr_parc;
            recon_mse_parc_version2(sbj,:)  = recon_mse_parc;
            recon_corr_parc_SC_exceptions_version2(sbj,:,:) = recon_corr_parc_SC_exceptions;
            recon_mse_parc_SC_exceptions_version2(sbj,:,:) = recon_mse_parc_SC_exceptions;

        elseif connectome_version == 3
            recon_corr_parc_version3(sbj,:) = recon_corr_parc;
            recon_mse_parc_version3(sbj,:)  = recon_mse_parc;
            recon_corr_parc_SC_exceptions_version3(sbj,:,:) = recon_corr_parc_SC_exceptions;
            recon_mse_parc_SC_exceptions_version3(sbj,:,:) = recon_mse_parc_SC_exceptions;

        elseif connectome_version == 4
            recon_corr_parc_version4(sbj,:) = recon_corr_parc;
            recon_mse_parc_version4(sbj,:)  = recon_mse_parc;
            recon_corr_parc_SC_exceptions_version4(sbj,:,:) = recon_corr_parc_SC_exceptions;
            recon_mse_parc_SC_exceptions_version4(sbj,:,:) = recon_mse_parc_SC_exceptions;
        end
    end
    
    if connectome_version == 1
        saveFileName = 'Results/fMRI_parc_reconstruction_all_subjects_EDRbinary';
        parsave2(saveFileName, recon_corr_parc_version1,recon_mse_parc_version1,...
        recon_corr_parc_SC_exceptions_version1,recon_mse_parc_SC_exceptions_version1)
    elseif connectome_version == 2
        saveFileName = 'project_laplacian/Results/fMRI_parc_reconstruction_all_subjects_EDRcontinuous';
        parsave2(saveFileName, recon_corr_parc_version2,recon_mse_parc_version2,...
        recon_corr_parc_SC_exceptions_version2,recon_mse_parc_SC_exceptions_version2)
    elseif connectome_version == 3
        saveFileName = 'project_laplacian/Results/fMRI_parc_reconstruction_all_subjects_Geometry';
        parsave2(saveFileName, recon_corr_parc_version3,recon_mse_parc_version3,...
        recon_corr_parc_SC_exceptions_version3,recon_mse_parc_SC_exceptions_version3)
    elseif connectome_version == 4
        saveFileName = 'project_laplacian/Results/fMRI_parc_reconstruction_all_subjects_EDRLR';
        parsave2(saveFileName, recon_corr_parc_version4,recon_mse_parc_version4,...
        recon_corr_parc_SC_exceptions_version4,recon_mse_parc_SC_exceptions_version4)
    end

end

elseif parLoop == 0

    %% load precalculated data
    
load('Results/fMRI_parc_reconstruction_all_subjects_EDRbinary.mat')

recon_corr_parc_version_all(1,:,:) = recon_corr_parc_version;
recon_mse_parc_version_all(1,:,:)  = recon_mse_parc_version;
recon_mse_parc_SC_exceptions_version_all(1,:,:,:)  = recon_mse_parc_SC_exceptions_version;
recon_corr_parc_SC_exceptions_version_all(1,:,:,:)  = recon_corr_parc_SC_exceptions_version;

load('Results/fMRI_parc_reconstruction_all_subjects_EDRcontinuous.mat')

recon_corr_parc_version_all(2,:,:) = recon_corr_parc_version;
recon_mse_parc_version_all(2,:,:)  = recon_mse_parc_version; 

recon_mse_parc_SC_exceptions_version_all(2,:,:,:)  = recon_mse_parc_SC_exceptions_version;
recon_corr_parc_SC_exceptions_version_all(2,:,:,:)  = recon_corr_parc_SC_exceptions_version;

load('Results/fMRI_parc_reconstruction_all_subjects_Geometry.mat')

recon_corr_parc_version_all(3,:,:) = recon_corr_parc_version;
recon_mse_parc_version_all(3,:,:)  = recon_mse_parc_version;

recon_mse_parc_SC_exceptions_version_all(3,:,:,:)  = recon_mse_parc_SC_exceptions_version;
recon_corr_parc_SC_exceptions_version_all(3,:,:,:)  = recon_corr_parc_SC_exceptions_version;

load('Results/fMRI_parc_reconstruction_all_subjects_EDRLR.mat')

recon_corr_parc_version_all(4,:,:) = recon_corr_parc_version;
recon_mse_parc_version_all(4,:,:)  = recon_mse_parc_version;

recon_mse_parc_SC_exceptions_version_all(4,:,:,:)  = recon_mse_parc_SC_exceptions_version;
recon_corr_parc_SC_exceptions_version_all(4,:,:,:)  = recon_corr_parc_SC_exceptions_version;
end
%% Figure 2A - fMRI correlation and MSE between reconstructed functional LR connections FC and empirical functional LR connections FC

figure
plot(squeeze(mean(recon_corr_parc_SC_exceptions_version_all(1,:,:,:),2))','Linewidth',2)
axis square;grid on;ylim([0, 1])
legend({'EDR binary','EDR weighted','Geometry','EDR+LR'},'Location','southeast')
ylabel('Correlation');xlabel('Modes')
%% Figure 2B
figure
violinplot(squeeze(recon_corr_parc_SC_exceptions_version_all(1,:,:,200)))
axis square;grid on;%ylim([0, 1])
set(gca,'XTick',1:4,'XTickLabel',{'EDR binary','EDR weighted','Geometry','EDR+LR'})
ylabel('Correlation')
%% Statistics for 2B
stats_modes_th = (squeeze(recon_corr_parc_SC_exceptions_version_all(1,:,:,200)));
[h_EDRbinary_EDRLR p_EDRbinary_EDRLR]               = ttest(stats_modes_th(:,1),stats_modes_th(:,4));
[h_EDRContinuous_EDRLR p_EDRcontinuous_EDRLR]       = ttest(stats_modes_th(:,2),stats_modes_th(:,4));
[h_Geometry_EDRLR p_Geometry_EDRLR]                 = ttest(stats_modes_th(:,3),stats_modes_th(:,4));
[h_Geometry_EDRContinuous p_Geometry_EDRContinuous] = ttest(stats_modes_th(:,2),stats_modes_th(:,3));

% %% different thesholds of functional LR connections (Correlation)
% figure
% for i=1:4
%     subplot(2,2,i)
%     plot(squeeze(mean(recon_corr_parc_SC_exceptions_version_all(i,:,:,:),2))','Linewidth',2)
%     title(strcat(['LR - euc.dist. > ',num2str(th_dis(i)),'mm and ', 'FC > 0.5']))
%     legend({'EDR binary','EDR weighted','Geometry','EDR+LR'},'Location','southeast')
%     axis square;grid on;ylim([0, 1]);ylabel('Correlation');xlabel('Modes')
% end
% %% fMRI correlation and MSE between reconstructed FC and empirical FC
% 
% figure,subplot(1,2,1);
% plot(squeeze(mean(recon_corr_parc_version_all,2))','Linewidth',2)
% legend({'EDR binary','EDR weighted','Geometry','EDR+LR'},'Location','southeast')
% ylabel('Correlation');xlabel('Modes');axis square
% subplot(1,2,2)
% plot(squeeze(mean(recon_mse_parc_version_all,2))','Linewidth',2)
% legend({'EDR binary','EDR weighted','Geometry','EDR+LR'},'Location','northeast')
% ylabel('MSE');xlabel('Modes');axis square
% 
% %% different thesholds of functional LR connections (MSE)
% figure
% for i=1:4
%     subplot(2,2,i)
%     plot(squeeze(mean(recon_mse_parc_SC_exceptions_version_all(i,:,:,:),2))','Linewidth',2)
%     title(strcat(['LR - euc.dist. > ',num2str(th_dis(i)),'mm and ', 'FC > 0.5']))
%     legend({'EDR binary','EDR weighted','Geometry','EDR+LR'},'Location','southeast')
%     axis square;grid on;ylabel('MSE');xlabel('Modes');ylim([0, 0.1])
% end
% 
% %% different thesholds of functional LR connections (Correlation) for mode = 200
% figure
% for i=1:4   
%     subplot(2,2,i)
%     violinplot(squeeze(recon_corr_parc_SC_exceptions_version_all(i,:,:,200)))
%     axis square;grid on;
%     set(gca,'XTick',1:4,'XTickLabel',{'EDR binary','EDR weighted','Geometry','EDR+LR'})
%     title(strcat(['LR - euc.dist. > ',num2str(th_dis(i)),'mm and ', 'FC > 0.5']))
%     ylabel('Correlation')
%
% end