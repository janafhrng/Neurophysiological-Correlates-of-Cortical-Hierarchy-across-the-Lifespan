%% PLOT FREQBANDS TOPOGRAPHIES AVERAGED
filepath = pwd;
path = [filepath '\PLOTS\'];

load([filepath '\DATA\03_specData_freqbands_norm']);

freq(1,:) = mean(zscore(delta_norm,[],2),1);
freq(2,:) = mean(zscore(theta_norm,[],2),1);
freq(3,:) = mean(zscore(alpha_norm,[],2),1);
freq(4,:) = mean(zscore(beta_norm,[],2),1);
freq(5,:) = mean(zscore(gamma_norm,[],2),1);
%% PLOT
parcellation = fetch_parcellation('fslr32k', 'schaefer', 200);
[surf_lh, surf_rh] = fetch_template_surface('fslr32k', 'layer', 'inflated');
load("vik.mat")

% delta
data2plot = freq(1,:)';
vertex_data = nan(size(parcellation));
for i = 1:200
    vertex_data(parcellation == i) = data2plot(i);
end

obj = plot_hemispheres(vertex_data, {surf_lh, surf_rh});
set(gcf,'units','points','position',[65,-65,1215,1215])
colormap([0.7 0.7 0.7;vik])
saveas(gcf,'x.png')
I1 = imread('x.png');
T1 = imcrop(I1, [180, 60, 555, 210]); % only left hemisphere
Imerged = imtile({T1},  'BorderSize', 4, 'BackgroundColor', 'white','GridSize', [1,1]);
imwrite(Imerged,[filepath '\PLOTS\topo_delta_uncorrected.png']);

% theta
data2plot = freq(2,:)';
vertex_data = nan(size(parcellation));
for i = 1:200
    vertex_data(parcellation == i) = data2plot(i);
end

obj = plot_hemispheres(vertex_data, {surf_lh, surf_rh});
set(gcf,'units','points','position',[65,-65,1215,1215])
colormap([0.7 0.7 0.7;vik])
saveas(gcf,'x.png')
I1 = imread('x.png');
T1 = imcrop(I1, [180, 60, 555, 210]); % only left hemisphere
Imerged = imtile({T1},  'BorderSize', 4, 'BackgroundColor', 'white','GridSize', [1,1]);
imwrite(Imerged,[filepath '\PLOTS\topo_theta_uncorrected.png']);

% alpha
data2plot = freq(3,:)';
vertex_data = nan(size(parcellation));
for i = 1:200
    vertex_data(parcellation == i) = data2plot(i);
end

obj = plot_hemispheres(vertex_data, {surf_lh, surf_rh});
set(gcf,'units','points','position',[65,-65,1215,1215])
colormap([0.7 0.7 0.7;vik])
saveas(gcf,'x.png')
I1 = imread('x.png');
T1 = imcrop(I1, [180, 60, 555, 210]); % only left hemisphere
Imerged = imtile({T1},  'BorderSize', 4, 'BackgroundColor', 'white','GridSize', [1,1]);
imwrite(Imerged,[filepath '\PLOTS\topo_alpha_uncorrected.png']);

% beta
data2plot = freq(4,:)';
vertex_data = nan(size(parcellation));
for i = 1:200
    vertex_data(parcellation == i) = data2plot(i);
end

obj = plot_hemispheres(vertex_data, {surf_lh, surf_rh});
set(gcf,'units','points','position',[65,-65,1215,1215])
colormap([0.7 0.7 0.7;vik])
saveas(gcf,'x.png')
I1 = imread('x.png');
T1 = imcrop(I1, [180, 60, 555, 210]); % only left hemisphere
Imerged = imtile({T1},  'BorderSize', 4, 'BackgroundColor', 'white','GridSize', [1,1]);
imwrite(Imerged,[filepath '\PLOTS\topo_beta_uncorrected.png']);

% gamma
data2plot = freq(5,:)';
vertex_data = nan(size(parcellation));
for i = 1:200
    vertex_data(parcellation == i) = data2plot(i);
end

obj = plot_hemispheres(vertex_data, {surf_lh, surf_rh});
set(gcf,'units','points','position',[65,-65,1215,1215])
colormap([0.7 0.7 0.7;vik])
saveas(gcf,'x.png')
I1 = imread('x.png');
T1 = imcrop(I1, [180, 60, 555, 210]); % only left hemisphere
Imerged = imtile({T1},  'BorderSize', 4, 'BackgroundColor', 'white','GridSize', [1,1]);
imwrite(Imerged,[filepath '\PLOTS\topo_gamma_uncorrected.png']);




%% PLOT TOPOGRAPHIES FOOOF EXP AND COG

% plot fooof exponent
% load([filepath '\DATA\foof_variables.mat'])
% fooof_plot = mean(zscore(fooof_exp_2_60Hz,[],2),1);
% we look at the fooof exponent from 2-60 Hz, the other one is from 1-60 Hz

% plot_hemispheres(fooof_plot', {surf_lh,surf_rh}, ...
%     'parcellation',labeling.schaefer_200)
% load('vik.mat')
% colormap([0.7 0.7	0.7;vik])
% set(gcf,'units','points','position',[50,-50,1200,1200])
% saveas(gcf,['x.png'])
% images = imcrop(imread(['x.png']), [170, 30, 1105, 270]);
% imwrite(images,[path, 'Tabl1_topo_FOOOFexp.png']);

% plot peak frrequency: center of gravity
load([filepath '\DATA\peak_freq_cog.mat'])
peakfreq_plot = mean(zscore(peak_freq_cog,[],2),1);

data2plot = peakfreq_plot';
vertex_data = nan(size(parcellation));
for i = 1:200
    vertex_data(parcellation == i) = data2plot(i);
end

obj = plot_hemispheres(vertex_data, {surf_lh, surf_rh});
set(gcf,'units','points','position',[65,-65,1215,1215])
colormap([0.7 0.7 0.7;vik])
saveas(gcf,'x.png')
I1 = imread('x.png');
T1 = imcrop(I1, [180, 60, 555, 210]); % only left hemisphere
Imerged = imtile({T1},  'BorderSize', 4, 'BackgroundColor', 'white','GridSize', [1,1]);
imwrite(Imerged,[filepath '\PLOTS\topo_cog.png']);


% %% PLOT StatAvl100
% load([filepath '\DATA\mat_z.mat'])
% 
% % create average over the StatAvl100 feature #448
% data2plot = squeeze(mean(mat_z(:,:,448),1));
% 
% plot_hemispheres(data2plot', {surf_lh,surf_rh}, ...
%     'parcellation',labeling.schaefer_200)
% load('vik.mat')
% colormap([0.7 0.7	0.7;vik])
% set(gcf,'units','points','position',[50,-50,1200,1200])
% saveas(gcf,['x.png'])
% images = imcrop(imread(['x.png']), [170, 30, 1105, 270]);
% imwrite(images,[path, 'topo_StatAvl100.png']);




