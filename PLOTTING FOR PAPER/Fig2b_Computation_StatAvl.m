% plot the data for a figure how to compute statavl100
filepath = pwd;

% plots for the features: DN_RemovePoints_absclose_08_remove_ac2diff and StatAvl100
load('lajolla10.mat')
load('vik.mat')

% adjust the data according to our schäfer 200 parcel structure
load([filepath '\DATA\mapping_idx_camcan2schaefer.mat'])

% loading the data of the respective participant --> cannot be done in the published code, for this access to the CAMCan repository is neccessary
% listing all participants files
Files = dir(fullfile([filepath '\raw_data\', '*singleshell.mat']));

% finding parcels which are low, medium or high in StatAvl100
load([filepath '\DATA\mat_z.mat'])
% loading the feature file with the complete list of features
load([filepath '\DATA\features_selected.mat']);
features= labels_selected;
% find the feature index
index = find(strcmp(features.Name,'StatAvl100'));
% we take sbj == 1 and find the parcels with the highest, medium and lowest
% value in StatAvl100
parcels2sort = squeeze(mat_z(1,:,index));

idx_h_min = find(parcels2sort == min(parcels2sort)); 
idx_h_med = find(parcels2sort == min(abs(parcels2sort))); 
idx_h_max = find(parcels2sort == max(parcels2sort));  
hierarchy_labels = {'idx_h_min','idx_h_med','idx_h_max'};


%% plotting the parcels
% welches sind überhaupt die parcels, die ich mir anschaue? 

%addpath(genpath('C:\Users\fehring\sciebo\Toolboxes'))
[surf_lh, surf_rh] = load_conte69();
labeling = load_parcellation('schaefer',200);

colors = [vik(20,:);vik(180,:);vik(230,:)];
data2plot = zeros(1,200);
for parcel = 1:3
    data2plot(eval(hierarchy_labels{parcel})) = parcel;
end
obj = plot_hemispheres(data2plot', {surf_lh,surf_rh},'parcellation',labeling.schaefer_200);
%set(gcf,'units','points','position',[50,-50,1200,1200])
colormap([0.7 0.7 0.7; colors])


%% NOT POSSIBLE IN THE PUBLISHED CODE
%% load one data file
subj = 2;
tmp_name = Files(subj).name;
load([filepath '\raw_data\' tmp_name]);
ts = ts(tmp_idx,:,:); % data: for each participant: parcels(214) x trials(30) x samples(3000) --> only use the schäfer 200 parcels

for i = 1:4; ticks{i} = num2str(round(i*100/300*1000));end

figure('color','w','units','normalized','outerposition',[0 0 1 1])
subplot(3,1,1)
data2plot = squeeze(ts(idx_h_min,2,1:400));
plot(data2plot,'-','LineWidth',2,'Color',vik(20,:))
hold on
ylim([-25 25])
set(gca,'YColor','none')
box off
for i = 1:4
    M = mean(data2plot((i-1)*100+1:i*100));
    plot([(i-1)*100+1, i*100],[M, M],'Color',[.5 .5 .5 .5],'LineWidth',2)
end
axis off
%set(gca,'FontSize',30)
hold off

subplot(3,1,2)
data2plot = squeeze(ts(idx_h_med,2,1:400));
plot(data2plot,'-','LineWidth',2,'Color',vik(180,:))
hold on
set(gca,'YColor','none')
ylim([-25 25])
box off
for i = 1:4
    M = mean(data2plot((i-1)*100+1:i*100));
    plot([(i-1)*100+1, i*100],[M, M],'Color',[.5 .5 .5 .5],'LineWidth',2)
end
axis off
%set(gca,'FontSize',30)
hold off

subplot(3,1,3)
data2plot = squeeze(ts(idx_h_max,2,1:400));
plot(data2plot,'-','LineWidth',2,'Color',vik(230,:))
hold on
set(gca,'YColor','none')
ylim([-25 25])
box off
for i = 1:4
    M = mean(data2plot((i-1)*100+1:i*100));
    plot([(i-1)*100+1, i*100],[M, M],'Color',[.5 .5 .5 .5],'LineWidth',2)
end
xticks([100 200 300 400])
xticklabels([ticks])
set(gca,'FontSize',30)
hold off

exportgraphics(gcf,[filepath '\PLOTS\Fig2b_computation_StatAvl_100.png'],'Resolution',500)













