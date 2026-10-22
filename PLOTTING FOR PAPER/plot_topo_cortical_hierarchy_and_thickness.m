%% PLOT CORTICAL HIERARCHY AND CORTICAL THICKNESS
filepath = pwd;

load('vik.mat')

%addpath(genpath('C:\Users\fehring\sciebo\Toolboxes'))
parcellation = fetch_parcellation('fsaverage5', 'schaefer', 200);
[surf_lh, surf_rh] = fetch_template_surface('fsaverage5', 'layer', 'inflated');


%% plotting hierarchy
load([filepath '\DATA\hierarchy_parc.csv'])

data2plot = hierarchy_parc';
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
imwrite(Imerged,[filepath '\PLOTS\topo_cortical_hierarchy.png']);


%% plotting cortical thickness

% loading cortical thickness
load([filepath '\DATA\thickness_smoothed_12_schaefer200_reordered.mat']);
thickness(306,:) = [];

thick_mean = mean(thickness,1);

data2plot = thick_mean';
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
imwrite(Imerged,[filepath '\PLOTS\topo_cortical_thickness.png']);
