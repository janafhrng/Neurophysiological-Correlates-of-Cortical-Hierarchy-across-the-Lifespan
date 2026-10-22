clearvars;
close all
clc

filepath = pwd;


%% plotting topogrpahies of fooof offset and exponent
% with knee!

load([filepath '\DATA\all_subjects_fooof_exp_knee.mat'])
load([filepath '\DATA\all_subjects_fooof_offs_knee.mat'])

% load index of parcels to map to schäfer200
load([filepath '\DATA\mapping_idx_camcan2schaefer']);


% EXPONENT --------------------------------------------------------------
% now we restucture this into 1 matrix containing: 
% parciticpants x parcels x features
exp  = zeros(350,214);
for i = 1:350 
    exp(i,:) = all_subjects_exponent{i};
end

% change the sorting of the data matrix from fsaverage to schäfer200
% look at the matrix only at the schäfer200 positions
exp = exp(:,tmp_idx);



% OFFSET ---------------------------------------------------------------
% now we restucture this into 1 matrix containing: 
% parciticpants x parcels x features
offs  = zeros(350,214);
for i = 1:350 
    offs(i,:) = all_subjects_offset{i};
end

% change the sorting of the data matrix from fsaverage to schäfer200
% look at the matrix only at the schäfer200 positions
offs = offs(:,tmp_idx);
save([filepath '\DATA\correctedFooofParam_fooofWknee.mat'],'exp',"offs");



%% z-score
exp_z = zscore(exp,[],2);
offs_z = zscore(offs,[],2);

save([filepath '\DATA\correctedFooofParam_fooofWknee_z.mat'],'exp_z',"offs_z");


%% plotting like before
% plot
parcellation = fetch_parcellation('fslr32k', 'schaefer', 200);
[surf_lh, surf_rh] = fetch_template_surface('fslr32k', 'layer', 'inflated');

data2plot = [mean(exp_z,1)];
vertex_data = nan(size(parcellation));
for i = 1:200
    vertex_data(parcellation == i) = data2plot(i);
end

obj = plot_hemispheres(vertex_data, {surf_lh, surf_rh});
set(gcf,'units','points','position',[65,-65,1215,1215])
colormap([0.7 0.7 0.7; vik])
saveas(gcf,'x.png')
I1 = imread('x.png');
T1 = imcrop(I1, [180, 60, 555, 210]); % only left hemisphere
Imerged = imtile({T1},  'BorderSize', 4, 'BackgroundColor', 'white','GridSize', [1,1]);
imwrite(Imerged,[filepath '\PLOTS\corrected_fooofexponent_fooofWknee.png']);



data2plot = [mean(offs_z,1)];
vertex_data = nan(size(parcellation));
for i = 1:200
    vertex_data(parcellation == i) = data2plot(i);
end

obj = plot_hemispheres(vertex_data, {surf_lh, surf_rh});
set(gcf,'units','points','position',[65,-65,1215,1215])
colormap([0.7 0.7 0.7; vik])
saveas(gcf,'x.png')
I1 = imread('x.png');
T1 = imcrop(I1, [180, 60, 555, 210]); % only left hemisphere
Imerged = imtile({T1},  'BorderSize', 4, 'BackgroundColor', 'white','GridSize', [1,1]);
imwrite(Imerged,[filepath '\PLOTS\corrected_fooofoffset_fooofWknee.png']);




%% LMEM with hierarchy, age --> are the t-values changed a lot due to 

% loading hierarchy
load([filepath '\DATA\hierarchy_parc.csv'])

% loading cortical thickness
load([filepath '\DATA\thickness_smoothed_12_schaefer200_reordered.mat']);
thickness(306,:) = []; 

% laoding age
demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));  

% EXPONENT -----------------------------------------------------
% HIERARCHY ------------------------
subj=[];
parcel=[];
feat=[];
hierarchy=[];
age=[];
for k1=1:350 %subjects
    subj=[subj ones(1,200)*k1];
    parcel = [parcel 1:200];
    feat=[feat exp_z(k1,:)];
    hierarchy=[hierarchy hierarchy_parc'];
    age = [age ones(1,200)*age_years(k1)];
end
varnames={'subj','feat','parcel','hierarchy','age'};
tbl=table(subj',feat',parcel',hierarchy',age','VariableNames',varnames);
model= 'feat ~ age*hierarchy + (1|subj) + (1|parcel)';

try
    lme = fitlme(tbl,model);
    allst = lme.Coefficients.tStat;
    pval = lme.Coefficients.pValue;
end


% THICKNESS ------------------------
exp_z(306,:)       = [];
age_years(306,:) = [];  

subj=[];
parcel=[];
feat=[];
thick=[];
age=[];
for k1=1:349 %subjects
    subj=[subj ones(1,200)*k1];
    parcel = [parcel 1:200];
    feat=[feat exp_z(k1,:)];
    thick=[thick thickness(k1,:)];
    age = [age ones(1,200)*age_years(k1)];
end
varnames={'subj','parcel','feat','thick','age'};
tbl=table(subj',parcel',feat',thick',age','VariableNames',varnames);
model= 'feat ~ age*thick + (1|subj) + (1|parcel)';

try
    lme = fitlme(tbl,model);
    allst=lme.Coefficients.tStat;
    pval=lme.Coefficients.pValue;
end




% OFFSET ----------------------------------------------------------------
% loading hierarchy
load([filepath '\DATA\hierarchy_parc.csv'])

% loading cortical thickness
load([filepath '\DATA\thickness_smoothed_12_schaefer200_reordered.mat']);
thickness(306,:) = []; 

% laoding age
demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));  

% HIERARCHY ------------------------
subj=[];
parcel=[];
feat=[];
hierarchy=[];
age=[];
for k1=1:350 %subjects
    subj=[subj ones(1,200)*k1];
    parcel = [parcel 1:200];
    feat=[feat offs_z(k1,:)];
    hierarchy=[hierarchy hierarchy_parc'];
    age = [age ones(1,200)*age_years(k1)];
end
varnames={'subj','parcel','feat','hierarchy','age'};
tbl=table(subj',parcel',feat',hierarchy',age','VariableNames',varnames);
model= 'feat ~ age*hierarchy + (1|subj) + (1|parcel)';

try
    lme = fitlme(tbl,model);
    allst = lme.Coefficients.tStat;
    pval = lme.Coefficients.pValue;
end


% THICKNESS ------------------------
offs_z(306,:)       = [];
age_years(306,:) = [];  

subj=[];
parcel=[];
feat=[];
thick=[];
age=[];
for k1=1:349 %subjects
    subj=[subj ones(1,200)*k1];
    parcel = [parcel 1:200];
    feat=[feat offs_z(k1,:)];
    thick=[thick thickness(k1,:)];
    age = [age ones(1,200)*age_years(k1)];
end
varnames={'subj','parcel','feat','thick','age'};
tbl=table(subj',parcel',feat',thick',age','VariableNames',varnames);
model= 'feat ~ age*thick + (1|subj) + (1|parcel)';

try
    lme = fitlme(tbl,model);
    allst=lme.Coefficients.tStat;
    pval=lme.Coefficients.pValue;
end














%% WITHOUT KNEE plotting topogrpahies of fooof offset and exponent
% without knee!

load([filepath '\DATA\all_subjects_fooof_exp_noknee_Feb26.mat'])
load([filepath '\DATA\all_subjects_fooof_offs_noknee_Feb26.mat'])

% load index of parcels to map to schäfer200
load([filepath '\DATA\mapping_idx_camcan2schaefer']);

% EXPONENT --------------------------------------------------------------
% now we restucture this into 1 matrix containing: 
% parciticpants x parcels x features
exp  = zeros(350,214);
for i = 1:350 
    exp(i,:) = all_subjects_exponent{i};
end

% change the sorting of the data matrix from fsaverage to schäfer200
% look at the matrix only at the schäfer200 positions
exp = exp(:,tmp_idx);

%exp(56:58,:) = NaN;

% plot
parcellation = fetch_parcellation('fslr32k', 'schaefer', 200);
[surf_lh, surf_rh] = fetch_template_surface('fslr32k', 'layer', 'inflated');

data2plot = mean(exp,1,'omitnan');
vertex_data = nan(size(parcellation));
for i = 1:200
    vertex_data(parcellation == i) = data2plot(i);
end

obj = plot_hemispheres(vertex_data, {surf_lh, surf_rh});




% OFFSET ---------------------------------------------------------------
% now we restucture this into 1 matrix containing: 
% parciticpants x parcels x features
offs  = zeros(350,214);
for i = 1:350 
    offs(i,:) = all_subjects_offset{i};
end

% change the sorting of the data matrix from fsaverage to schäfer200
% look at the matrix only at the schäfer200 positions
offs = offs(:,tmp_idx);

% plot
parcellation = fetch_parcellation('fslr32k', 'schaefer', 200);
[surf_lh, surf_rh] = fetch_template_surface('fslr32k', 'layer', 'inflated');

data2plot = mean(offs,1);
vertex_data = nan(size(parcellation));
for i = 1:200
    vertex_data(parcellation == i) = data2plot(i);
end

obj = plot_hemispheres(vertex_data, {surf_lh, surf_rh});













