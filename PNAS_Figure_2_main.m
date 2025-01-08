%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 2
%%% Original: James Pang, Monash University, 
%%% updated: Jakub Vohryzek, Universitat Pompeu Fabra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% loading data
% change directory here
repo_dir = '/BrainEigenmodes_EDRLR-main';
addpath(genpath(repo_dir))

%% Load surface files for visualization
hemisphere = 'lh';
num_modes = 200;
num_sbj = 255 %
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

th_dis = 40;

%%

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
    %% pre-allocation

    recon_mse_parc_version1 = zeros(num_sbj, num_modes);
    recon_mse_parc_version2 = zeros(num_sbj, num_modes);
    recon_mse_parc_version3 = zeros(num_sbj, num_modes);
    recon_mse_parc_version4 = zeros(num_sbj, num_modes);

    recon_mse_parc_SC_exceptions_version1 = zeros(num_sbj,size(th_dis,2), num_modes);
    recon_mse_parc_SC_exceptions_version2 = zeros(num_sbj,size(th_dis,2), num_modes);
    recon_mse_parc_SC_exceptions_version3 = zeros(num_sbj,size(th_dis,2), num_modes);
    recon_mse_parc_SC_exceptions_version4 = zeros(num_sbj,size(th_dis,2), num_modes);
    %%
    for sbj = 1:num_sbj
        disp(['connectome v',num2str(i),' subject #',num2str(sbj)])
    
        % load all-subjects= rfMIR timeseries data
        data = load(['/Datasets/HCP/HCP255/subject_',num2str(subjectList(sbj)),sprintf('_rfMRI_timeseries-%s.mat', hemisphere)]);
        data_to_reconstruct = data.timeseries;
        T = size(data_to_reconstruct, 2);
        recon_beta = zeros(num_modes, T, num_modes);
        for mode = 1:num_modes
            basis = eigenmodes(cortex_ind, 1:mode);     
            recon_beta(1:mode,:,mode) = calc_eigendecomposition(data_to_reconstruct(cortex_ind,:), basis, 'matrix');
        end

        % Calculate empirical FC
        data_parc_emp = calc_parcellate(parc, data_to_reconstruct);
        data_parc_emp = calc_normalize_timeseries(data_parc_emp');
        data_parc_emp(isnan(data_parc_emp)) = 0;
        
        FC_emp_temp = data_parc_emp'*data_parc_emp;
        FC_emp_temp = FC_emp_temp/T; % empirical FC
        FCvec_emp = FC_emp_temp(triu_ind); % this is what is being fitted
    
        % Calculate reconstructed FC and accuracy (slow to run with more modes)
        FCvec_recon = zeros(length(triu_ind), num_modes);
        recon_mse_parc = zeros(1, num_modes);               
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
                   
            recon_mse_parc(mode) = immse(FCvec_emp, FCvec_recon(:,mode));

            % functional long-range exceptions
            for t = 1:4
                FC_LRE_matrix = ((FC_LRE_euc_dist_matrix>th_dis(t)).*(FC_emp_HCP_255_subject_mean>0.5)); 
                a = FC_LRE_matrix.*FC_emp_temp;
                b = FC_LRE_matrix.*FC_recon_temp;
                recon_mse_parc_SC_exceptions(t,mode)   = immse(nonzeros(a(triu_ind)), nonzeros(b(triu_ind)));

            end
        end
    
        if connectome_version == 1
            recon_mse_parc_version1(sbj,:)  = recon_mse_parc;        
            recon_mse_parc_SC_exceptions_version1(sbj,:,:) = recon_mse_parc_SC_exceptions;

        elseif connectome_version == 2
            recon_mse_parc_version2(sbj,:)  = recon_mse_parc;
            recon_mse_parc_SC_exceptions_version2(sbj,:,:) = recon_mse_parc_SC_exceptions;

        elseif connectome_version == 3
            recon_mse_parc_version3(sbj,:)  = recon_mse_parc;
            recon_mse_parc_SC_exceptions_version3(sbj,:,:) = recon_mse_parc_SC_exceptions;

        elseif connectome_version == 4
            recon_mse_parc_version4(sbj,:)  = recon_mse_parc;
            recon_mse_parc_SC_exceptions_version4(sbj,:,:) = recon_mse_parc_SC_exceptions;
        end
    end
    
    if connectome_version == 1
        saveFileName = 'Results/fMRI_parc_reconstruction_all_subjects_EDRbinary';
        parsave2(saveFileName,recon_mse_parc_version1,recon_mse_parc_SC_exceptions_version1)
    elseif connectome_version == 2
        saveFileName = 'project_laplacian/Results/fMRI_parc_reconstruction_all_subjects_EDRcontinuous';
        parsave2(saveFileName, recon_mse_parc_version2,recon_mse_parc_SC_exceptions_version2)
    elseif connectome_version == 3
        saveFileName = 'project_laplacian/Results/fMRI_parc_reconstruction_all_subjects_Geometry';
        parsave2(saveFileName,recon_mse_parc_version3,recon_mse_parc_SC_exceptions_version3)
    elseif connectome_version == 4
        saveFileName = 'project_laplacian/Results/fMRI_parc_reconstruction_all_subjects_EDRLR';
        parsave2(saveFileName,recon_mse_parc_version4,recon_mse_parc_SC_exceptions_version4)
    end

end

elseif parLoop == 0

%% load precalculated data
    
load('Results/long_3T/fMRI_parc_reconstruction_all_subjects_EDRbinary.mat')

recon_mse_parc_version_all(1,:,:)  = recon_mse_parc_version;
recon_mse_parc_SC_exceptions_version_all(1,:,:,:)  = recon_mse_parc_SC_exceptions_version;

load('Results/long_3T/fMRI_parc_reconstruction_all_subjects_EDRcontinuous.mat')

recon_mse_parc_version_all(2,:,:)  = recon_mse_parc_version; 

recon_mse_parc_SC_exceptions_version_all(2,:,:,:)  = recon_mse_parc_SC_exceptions_version;

load('Results/long_3T/fMRI_parc_reconstruction_all_subjects_Geometry.mat')

recon_mse_parc_version_all(3,:,:)  = recon_mse_parc_version;

recon_mse_parc_SC_exceptions_version_all(3,:,:,:)  = recon_mse_parc_SC_exceptions_version;

load('Results/long_3T/fMRI_parc_reconstruction_all_subjects_EDRLR.mat')

recon_mse_parc_version_all(4,:,:)  = recon_mse_parc_version;

recon_mse_parc_SC_exceptions_version_all(4,:,:,:)  = recon_mse_parc_SC_exceptions_version;
end
%% data curation
data2mse =  [[squeeze(recon_mse_parc_SC_exceptions_version_all(1,:,1,200)),...
    squeeze(recon_mse_parc_SC_exceptions_version_all(2,:,1,200)),...
    squeeze(recon_mse_parc_SC_exceptions_version_all(3,:,1,200)),...
    squeeze(recon_mse_parc_SC_exceptions_version_all(4,:,1,200))]'];
%% auxillary

group_names = {'3T Rest'};
condition_names= {'EDR binary', 'EDR continuous' , 'Geometry', 'EDR+LR'};

% an alternative color scheme for some plots
c =  [0.45, 0.80, 0.69;...
      0.98, 0.40, 0.35;...
      0.55, 0.60, 0.79;...
      0.90, 0.70, 0.30]; 
group_inx = [ones(1,255), 2.*ones(1,255) 3.*ones(1,255) 4.*ones(1,255)];


%% MSE

figure
% different color scheme, different position of boxplots and scatter
h = daviolinplot(data2mse(:,1:1),'groups',group_inx,'color',c,'outsymbol','k+',...
    'xtlabels', condition_names(1:4),'scatter',2,'jitter',1,...
    'box',1,'boxcolors','same','scattercolors','same',...
    'boxspacing',1.1);%,'legend',group_names(1));
ylabel('MSE');
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]); % make more space for the legend
set(gca,'FontSize',10); grid on

%% MSE
figure
plot(squeeze(mean(recon_mse_parc_SC_exceptions_version_all(:,:,1,:),2))','Linewidth',2)
axis square;grid on;ylim([0, 0.05])
legend({'EDR binary','EDR weighted','Geometry','EDR+LR'},'Location','northeast')
ylabel('MSE');xlabel('Modes');
set(gcf,'Color', [1 1 1]);

%% Statistics
stats_modes_th = (squeeze(recon_mse_parc_SC_exceptions_version_all(:,:,1,200)));
[h_EDRbinary_EDRLR, p_EDRbinary_EDRLR]               = ttest(stats_modes_th(1,:),stats_modes_th(4,:));
[h_EDRContinuous_EDRLR, p_EDRcontinuous_EDRLR]       = ttest(stats_modes_th(2,:),stats_modes_th(4,:));
[h_Geometry_EDRLR, p_Geometry_EDRLR]                 = ttest(stats_modes_th(3,:),stats_modes_th(4,:));
[h_Geometry_EDRContinuous, p_Geometry_EDRContinuous] = ttest(stats_modes_th(2,:),stats_modes_th(3,:));
P_corrected = ([p_EDRbinary_EDRLR,p_EDRcontinuous_EDRLR,p_Geometry_EDRLR,p_Geometry_EDRContinuous]).*4;
