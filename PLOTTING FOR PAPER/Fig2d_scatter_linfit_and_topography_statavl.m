filepath = pwd;

load([filepath '\DATA\hierarchy_parc.csv'])
load([filepath '\DATA\mat_z.mat'])
features = readtable([filepath '\DATA\features_z.csv']);

load('vik.mat')

%addpath(genpath('C:\Users\fehring\sciebo\Toolboxes'))
[surf_lh, surf_rh] = load_conte69();
labeling = load_parcellation('schaefer',200);

%% PLOTTING LINFIT FOR STATAVL100
% find the feature index
index = find(strcmp(features.Var2,'StatAvl100'));

% creating a linear fit
mat_avg_z = mean(squeeze(mat_z(:,:,index)),1);
lm_h = fitlm(hierarchy_parc,mat_avg_z);

% scatter plot
figure
plot(lm_h)
xlabel('cortical hierarchy')
ylabel('feature value')
Lines = findobj(gca, 'Type','Line');
CI = findobj('Type','line','DisplayName','Confidence bounds');
CIX = CI.XData;
CIY = CI.YData;
lm_intercept = table2array(lm_h.Coefficients(1,1));
lm_slope = table2array(lm_h.Coefficients(2,1));

% in a nicer plot
dist = lm_intercept+lm_slope*CIX - CIY;
x = [CIX CIX(end:-1:1)];
y = [CIY CIY(end:-1:1)+dist*2];

% sort hierarchy_parc according to the hierarchy
[~,idx] = sort(mat_avg_z,'descend');
hierarchy_parc_sort = hierarchy_parc(idx);
mat_avg_z_sort = mat_avg_z(idx);

% find color vector indices

f = lm_intercept+lm_slope*hierarchy_parc;
figure('Color','w','Position',[1,1,400,600])
scatter(hierarchy_parc_sort,mat_avg_z_sort,90,vik(40,:),'filled','o','MarkerFaceAlpha',.5)
hold on
plot(hierarchy_parc,f,'LineWidth',3,'Color',[0.5 0.5 0.5])
fill(x,y,[0.5 0.5 0.5],'FaceAlpha',.4,'EdgeColor','none')
xlabel('cortical hierarchy','FontSize',16)
ylabel('stationarity of the mean','FontSize',16)
set(gca,'box','off') 
set(gca,'FontSize',20)
ylim([-1.5 2.5])
exportgraphics(gcf,[filepath '\PLOTS\Fig2d_scatter_linfit_statAvl100.png'],'Resolution',500)

%% stats
lm_h



%% topo StatAvl100 
% Load the atlas from the CamCan study, schäfer200
[surf_lh, surf_rh] = load_conte69();
labeling = load_parcellation('schaefer',200);


obj = plot_hemispheres(mat_avg_z', {surf_lh,surf_rh}, ...
    'parcellation',labeling.schaefer_200)
load('vik.mat')
colormap([0.7 0.7	0.7;vik])
set(gcf,'units','points','position',[50,-50,1200,1200])
saveas(gcf,['x.png'])
images = imcrop(imread(['x.png']),[170, 30, 1105, 270]);
imwrite(images,[filepath '\PLOTS\topo_StatAvl100.png']);
 
