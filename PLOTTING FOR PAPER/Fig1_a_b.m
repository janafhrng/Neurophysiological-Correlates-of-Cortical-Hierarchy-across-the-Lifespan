%% plots brain gray

filepath = pwd;

% plot
parcellation = fetch_parcellation('fslr32k', 'schaefer', 200);
[surf_lh, surf_rh] = fetch_template_surface('fslr32k', 'layer', 'inflated');

data2plot = repmat(1,200,1);

% obj = plot_hemispheres(data2plot, {surf_lh,surf_rh}, ...
%     'parcellation',labeling.schaefer_200, 'views', 'lm');

vertex_data = nan(size(parcellation));
for i = 1:200
    vertex_data(parcellation == i) = data2plot(i);
end

obj = plot_hemispheres(vertex_data, {surf_lh, surf_rh});
set(gcf,'units','points','position',[65,-65,1215,1215])
colormap([0.9 0.9 0.9])
saveas(gcf,'x.png')
I1 = imread('x.png');
T1 = imcrop(I1, [180, 60, 555, 210]); % only left hemisphere
Imerged = imtile({T1},  'BorderSize', 4, 'BackgroundColor', 'white','GridSize', [1,1]);
imwrite(Imerged,[filepath '\PLOTS\Fig1a_brain_gray.png']);


%% plotting the time series
% load one data file
subj = 1;

% listing all participants files --> not possible in the published code, as we do not publish the raw data
Files = dir(fullfile([filepath '\raw_data\', '*singleshell.mat']));
tmp_name = Files(subj).name;
load([filepath '\raw_data\' tmp_name]);

% adjust the data according to our schäfer 200 parcel structure
load([filepath '\DATA\mapping_idx_camcan2schaefer.mat'])
ts = ts(tmp_idx,:,:); % data: for each participant: parcels(214) x trials(30) x samples(3000) --> only use the schäfer 200 parcels


figure('color','w','units','normalized','outerposition',[0 0 1 1])
subplot(3,1,1)
data2plot = squeeze(ts(1,1,1:1000));
plot(data2plot,'-','LineWidth',2,'Color',[.3 .3 .3])
hold on
ylim([-25 25])
set(gca,'YColor','none')
set(gca,'XColor','none')
box off
set(gca,'FontSize',30)
set(gca,'xtick',[])
set(gca,'ytick',[])
hold off
exportgraphics(gca,[filepath '\PLOTS\Fig1a_timeseries_grey.png'],'Resolution',500)


%% plotting the power spectrum

% read the power spectrum
load([filepath '\DATA\01_specData_all.mat']);
specData_all = specData_all(:,tmp_idx,:);

% take one subject and one aprcel
data2plot = squeeze(specData_all(1,1,:));

% plot
figure('Color','white','Position',[1,1,600,400])
plot(data2plot(2:70),'Color',[.5 .5 .5],'LineWidth',2.5)
%xlim([2, 60])
box off
set(gca,'xtick',[])
set(gca,'ytick',[])
exportgraphics(gcf,[filepath '\PLOTS\Fig1b_powerspectrum.png'],'Resolution',500)

%% plotting time-series features
%loading matrix with the data in the form of participant x parcels x features
load([filepath '\DATA\mat_z.mat'])

% take one participant and one parcel
data2plot = squeeze(mat_z(2,:,100:400));

load('lajolla.mat')

% plot as a heatplot
figure('Color','white','Position',[1,1,1000,400])
h = heatmap(data2plot,Colormap=flip(lajolla));
grid off
%h.ColorScaling = 'scaledcolumns';

exportgraphics(h,[filepath '\PLOTS\Fig1b_heatmap_hctsa_feat.png'],'Resolution',500)


