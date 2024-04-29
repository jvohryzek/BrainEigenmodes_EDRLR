%% loading data
% change directory here
repo_dir = '/Users/jakub/Matlab/Collaboration_Deco/project_laplacian/BrainEigenmodes_EDRLR-main';
addpath(genpath(repo_dir))

%%
% load structural connectome
load('empirical/S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat')
connectome = avgSC_L;clear avgSC_L
% load euclidean distance matrix
load('empirical/synthetic_distance_connectome','surface_dist')

%% Initialisation
linfunc = @(A, x)(A(1)*x+A(2));
options=optimset('MaxFunEvals',10000,'MaxIter',1000,'Display','off');

NPARCELLS = 29696; % number of vertices
NR = 400; % number of bins
NRini = 20; % min bins
NRfin = 380; % max bins

NSTD = 3; % number of standard deviations
DistRange = 40;

Isubdiag = find(tril(ones(NPARCELLS),-1)); % indices for the lower triangular matrix

% normalise connectome between 0 and 1
C = connectome;
C = C/max(max(C));
clear connectome

rr = surface_dist; % defining the euclidean distance between vertices
clear surface_dist
range = max(max(rr)); % max distance
delta = range/NR; % setting up the bin width

for i=1:NR
    xrange(i) = delta/2+delta*(i-1); % setting up the bin edges
end

%% SC

%% %% option 1
numind3 = zeros(1,NR);
sc_3 = zeros(1,NR);
sc_density = cell(1,NR);
sc_density_i = cell(1,NR);
sc_density_j = cell(1,NR);

index = floor(rr/delta)+1; % what bin a particular distance falls into
index(index > NR) = NR; % Handle the case when index is NR+1

for n =1:400
    n
    idx = find(index == n); % for indexing matric entries all together
    [idx_i,idx_j] = find(index == n); % for saving i and j

    sc_density{n} = C(idx);
    sc_density_i{n} = idx_i;
    sc_density_j{n} = idx_j;
    sc_3(n) = sum(C(idx));
    % for with zeros
    numind3(n)=size(idx,1);
    clear idx idx_i idx_j
end

sctotal=sc_3./numind3;
for k=1:NR
    xcoor(k)=xrange(k);
    ycoor(k)=log(sctotal(k));
    ycoor2(k)=sctotal(k);
end
%% linear fitting between for 400 bins
expfunc = @(A, x)(A(1)*exp(-A(2)*x));
options = optimset('MaxFunEvals',10000,'MaxIter',1000,'Display','off');
A0 = [0.15 0.18];
Afit = lsqcurvefit(expfunc,A0,xcoor(25:end),ycoor2(25:end),[-100 -100],[100 100],options);
yl = Afit(1)*exp(-Afit(2)*xcoor);
lambda = Afit(2)
Afit

%% choosing the long-range connections above 3 standard deviations
Clong=zeros(NPARCELLS,NPARCELLS);
Clongall=zeros(NPARCELLS,NPARCELLS);
Clongrall=zeros(NPARCELLS,NPARCELLS);
Clongr=zeros(NPARCELLS,NPARCELLS);
numexcSC=zeros(NR,1);
valexcSC=zeros(NR,1);

for i = NRini:NRfin
    i
    mv=nanmean(sc_density{i}); % mean of a bin
    st=nanstd(sc_density{i}); % std of a bin
    numexcSC(i)=length(find(sc_density{i}>mv+NSTD*st))/length(sc_density{i}); % ratio of number of conn above 3std
    valexcSC(i)=nanmean(sc_density{i}(find(sc_density{i}>mv+NSTD*st)))/nanmean(sc_density{i});% ratio of conn values above 3std
    
    for j=1:length(sc_density{i})

       Clongall(sc_density_i{i}(j),sc_density_j{i}(j)) = sc_density{i}(j);
       Clongrall(sc_density_i{i}(j),sc_density_j{i}(j)) = xrange(i);
    end
    ind_exc = find(sc_density{i}>mv+NSTD*st);
    for n=1:size(ind_exc,1)
        ind = ind_exc(n)
        if rr(sc_density_i{i}(ind),sc_density_j{i}(ind))>DistRange
            Clong(sc_density_i{i}(ind),sc_density_j{i}(ind))=sc_density{i}(ind);
            Clongr(sc_density_i{i}(ind),sc_density_j{i}(ind))=xrange(i);
        end
    end
   clear ind
end
%% normalising
Clongall_v0 = Clongall;
Clong_v0 = Clong;
C_v0 = C;

Clongall=Clongall/Afit(1);
Clong=Clong/Afit(1);
SC=C/Afit(1);
%%
for i=1:400
    val(i)=nanmean(sc_density{i});
    valsd(i)=nanstd(sc_density{i});
end

%% Creating the EDR+LR connectome
EDR_conn = Afit(1)*exp(-Afit(2)*rr);
EDR_LRE_tmp = EDR_conn;
EDR_LRE_tmp_idx = find(Clong_v0 > 0);
EDR_LRE_tmp(EDR_LRE_tmp_idx) = 0;
EDR_LRE = EDR_LRE_tmp+Clong_v0;
% EDR_LRE_v0 = EDR_conn+Clong_v0;
%% Plotting the Connectome, EDR
figure,
subplot(121);imagesc(C,[0 0.06]);axis square;title('connectome');colorbar
subplot(122);imagesc(EDR_conn,[0 0.06]);axis square;title('EDR');colorbar
colormap(flipud(bone))
%%
% plot figure with 
figure
shadedErrorBar(xcoor,val,valsd)
hold on
plot(xcoor,yl,'r')
ylabel('Fibre density');xlabel('Distance (mm)')
%% Plotting just the EDRLR connections
figure,
hold on
shadedErrorBar(xcoor,val,valsd)
hold on
plot(xcoor,yl,'r')
ylabel('Fibre density');xlabel('Distance (mm)')
plot(Clongrall(1:100:size(Clongrall,1),1:100:size(Clongrall,1)),Clongall_v0(1:100:size(Clongrall,1),1:100:size(Clongrall,1)),'.y')
hold on
plot(Clongr(1:100:size(Clongr,1),1:100:size(Clongr,1)),Clong_v0(1:100:size(Clongr,1),1:100:size(Clongr,1)),'.r')


%% Plotting just the EDR+LR connections
figure,
hold on
shadedErrorBar(xcoor,val,valsd)
hold on
plot(xcoor,yl,'r')
ylabel('Fibre density');xlabel('Distance (mm)')
plot(Clongrall(1:100:size(Clongrall,1),1:100:size(Clongrall,1)),EDR_LRE(1:100:size(EDR_conn,1),1:100:size(EDR_conn,1)),'.r')
xlim([1 180])
legend('SC','EDR','LR','EDR','LR')
%% Comparions EDR, EDR+LR and the Connectome
figure,
subplot(235);imagesc(Clong>0);axis square;title('Clong>3std>40mm');colorbar
subplot(231);imagesc(EDR_conn,[0 0.09]);axis square;title('EDR');colorbar
subplot(232);imagesc(EDR_LRE,[0 0.09]);axis square;title('EDR+LR');colorbar
subplot(233);imagesc(C,[0 0.09]);axis square;title('Connectome');colorbar
colormap(flipud(bone))
%% Comparions EDR+LR and Euclidean Distance
figure,
subplot(221);imagesc(Clong>0);axis square;title('Clong>3std>40mm')
subplot(222);imagesc(Clongall>0);axis square;title('C')
subplot(223);imagesc(Clongr);axis square;title('euc. dist. of Clong>3std>40mm')
subplot(224);imagesc(Clongrall);axis square;title('euc. dist.')
colormap(flipud(bone))
%% Downsampling the matrices for histogram plotting
C_donwsampled = C(1:100:size(Clongrall,1),1:100:size(Clongrall,1));
EDR_donwsampled = EDR_conn(1:100:size(Clongrall,1),1:100:size(Clongrall,1));
EDR_LR_donwsampled = EDR_LR(1:100:size(Clongrall,1),1:100:size(Clongrall,1));

%% Plotting the Connectome, EDR, EDR+LR histograms
figure,
subplot(131);hist(nonzeros(C_donwsampled),100);title('connectome');ylim([0 200]);xlim([0 0.06]);
subplot(132);hist(nonzeros(EDR_donwsampled),100);title('EDR');ylim([0 200]);xlim([0 0.06]);
subplot(133);hist(nonzeros(EDR_LR_donwsampled),100);title('EDR+LR');ylim([0 200]);xlim([0 0.06]);
