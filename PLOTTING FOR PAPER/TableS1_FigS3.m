%% LIST OF THE 50 TOP FEATURES IN LMEM FOR MAIN CORTICAL HIERARCHY
% loading the feature file with the complete list of features 
filepath = pwd;

features = readtable([filepath '\DATA\features_z.csv']);

% loading lme stats
allst = readmatrix([filepath '\DATA\lme_out_age_effect_hierarchy.csv']);

% ranking the features
% we rank the features according to the highes t-score in order to find
% features that show a really strong gradient along cortical hierarchy
t = allst(:,3);
[~,rank] = sort(abs(t),'descend');
t_ranked = t(rank);
features_ranked = features(rank,:);

data2plot = [features_ranked(1:50,"Var2"),array2table(t_ranked(1:50))];
writetable(data2plot,[filepath '\PLOTS\SUPPL_table_lmem_tValues_list.csv'])

%% LIST OF ALL STATAVL FEATURES AND THEIR T-VALUES
% feature index = 446:253

idx = 446:453;
data2plot = [features(idx,"Var2"),array2table(t(idx))];

% and plot
%addpath(genpath('C:\Users\fehring\sciebo\Toolboxes'))
[surf_lh, surf_rh] = load_conte69();
labeling = load_parcellation('schaefer',200);
load("vik.mat")

load([filepath '\DATA\mat_z.mat'])
data2plot = squeeze(mean(zscore(squeeze(mat_z(:,:,idx)),[],2),1));

obj = plot_hemispheres(data2plot(:,1:4), {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_200, 'views', 'lm');
set(gcf,'units','points','position',[50,-50,1200,1200])
colormap([0.7 0.7	0.7;vik])
saveas(gcf,'x.png')
I1 = imread('x.png');
T1 = imcrop(I1, [145, 30, 583, 240]); % only left hemisphere
T2 = imcrop(I1, [145, 315, 583, 240]);
T3 = imcrop(I1, [145, 600, 583, 240]);
T4 = imcrop(I1, [145, 885, 583, 240]);

Imerged = imtile({T1,T2,T3,T4},  'BorderSize', 4, 'BackgroundColor', 'white','GridSize', [4,1]);
imwrite(Imerged,[filepath '\PLOTS\SUPPL_table_StatAvl_windowSize_topo1.png']);


obj = plot_hemispheres(data2plot(:,5:8), {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_200, 'views', 'lm');
set(gcf,'units','points','position',[50,-50,1200,1200])
colormap([0.7 0.7	0.7;vik])
saveas(gcf,'x.png')
I1 = imread('x.png');
T1 = imcrop(I1, [145, 30, 583, 240]); % only left hemisphere
T2 = imcrop(I1, [145, 315, 583, 240]);
T3 = imcrop(I1, [145, 600, 583, 240]);
T4 = imcrop(I1, [145, 885, 583, 240]);
Imerged = imtile({T1,T2,T3,T4},  'BorderSize', 4, 'BackgroundColor', 'white','GridSize', [4,1]);
imwrite(Imerged,[filepath '\PLOTS\SUPPL_table_StatAvl_windowSize_topo1.png']);

