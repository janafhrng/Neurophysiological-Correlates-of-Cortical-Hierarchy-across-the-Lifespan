clearvars;
close all
clc

%% Here, we use the data, which has been already been preprocessed by Christina Stier
% She averaged for each participant over the all time segements, so we dont
% not have to this here anymore. Further, she already filtered out some
% data. However, to use the plotting functions later on, we have to reorder
% our data set and look at only the 200 schäfer200 parcels, disregarding
% the 14 subcortical parcels included here.

%% online or averaged coordinates
% it is possible to use the averaged coordinates across the vertizes of
% each parcel OR the coordinates that can be found online for each parcel.
% They are quite similar. 
% Importantly, the data and the averaged parcel coordinates are not structured as the
% schäfer 200 parcel coordinates. The rows, regarding the actual parcels
% are sorted differently. This will create problems when plotting the data
% topographically. Thus, we have to reorder
% - data: the data matrix has to be reordered for the averaged and online
% parcellation according to the way the schäfer 200 parcels are sorted.
% - parcel coordinates: the parcel coordinates only have to be reordered
% when we use the averaged parcel coordinates.

filepath = pwd;
 
%% load coordinates
% online coordinates
labels_online = readtable([filepath '\DATA\Schaefer2018_200Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv']);


% fsaverage coordinates to average across vertex
load([filepath '\DATA\Schaefer2018_200Parcels_7Networks_suma-all-fsaverage-10.mat'])

%% create averaged coordinates for parcels
% map vertices to the parcels we are interested in and average across
% vertices per parcel + save the name (key) of each parcel
label_pos = zeros(214,3);
label_key = {};
for i=1:214
    tmp_key = suma_all.annot_key{1,1}(i);
    tmp_roiidx=find(suma_all.annot==tmp_key);
    tmp_roiidx_pos=suma_all.pos(tmp_roiidx,:);
    tmp_pos = mean(tmp_roiidx_pos,1);
    label_pos(i,:) = tmp_pos;
    tmp_label_key = suma_all.annot_key{1,2}{i};
    label_key = [label_key; {tmp_label_key(1:end-3)}];
end

%% map the sorting of the fsaverage parcellation to the schäfer (online) parcellation
% find which row of the labels from fsaverage belong to which row of the
% schäfer atlas
for i = 1:200
    tmp_idx(i) = find(strcmp(label_key,labels_online.ROIName(i)));
end

% save tmp_idx
save([filepath '\DATA\mapping_idx_camcan2schaefer'],'tmp_idx');

%% load data matrix and average across time segments
% loading all data: cell structure with 1 cell per participant (350)
% with a matrix containing: parcels(214) x hcsta features (5961)
load([filepath '\DATA\all_subj_feat.mat']);

% now we restucture this into 1 matrix containing: 
% parciticpants x parcels x features
mat  = zeros(350,214,5961);
for i = 1:350 
    mat(i,:,:) = all_subj_feat{1,i};
end

% change the sorting of the data matrix from fsaverage to schäfer200
% look at the matrix only at the schäfer200 positions
mat = mat(:,tmp_idx,:);
save([filepath '\DATA\all_subj_feat_reordered_schafer200'], "mat", '-v7.3');

%% sort the whole power spectrum accordingly
load([filepath '\DATA\01_specData_all.mat'],'specData_all')
specData = specData_all(:,tmp_idx,:);
save([filepath '\DATA\02_specData_all_reordered_schaefer200.mat'], "specData");

%% sort cortical thinkness accordingly
load([filepath '\DATA\thickness_smoothed_12_schaefer200.mat'])
% reoder this matrix, so that in the end we have: participant(350) x parcel(214)
thickness  = zeros(350,214);

for i = 1:350 
    thickness(i,:) = schaefer_all{i,1};
end

thickness = thickness(:,tmp_idx);
save([filepath 'Z\DATA\thickness_smoothed_12_schaefer200_reordered'], "thickness");

%% change the sorting of centre of gravity (alpha peak frequency)
load([filepath '\DATA\all_subjects_paf_cog.mat'])

peak_freq_cog  = zeros(350,214);
for i = 1:350 
    peak_freq_cog(i,:) = all_subjects_cog{i,1};
end

peak_freq_cog = peak_freq_cog(:,tmp_idx);
writematrix(peak_freq_cog,[filepath '\DATA\peak_freq_cog.csv'])

%% change the sorting of FOOF variables (offset and exponent)
%The data contains on the FOOF offset and exponent for the fe data from
%1-60Hz AND 2-60Hz (these are marked differently)
load([filepath '\DATA\all_subjects_fooof_exp_2_60Hz.mat']);
all_subjects_exponent_2_60Hz = all_subjects_exponent;
load([filepath '\DATA\all_subjects_fooof_exp.mat'],'all_subjects_exponent');
load([filepath '\DATA\all_subjects_fooof_offs_2_60Hz.mat']);
all_subjects_offset_2_60Hz = all_subjects_offset;
load([filepath '\DATA\all_subjects_fooof_offs.mat']);

fooof_exp  = zeros(350,214);
fooof_exp_2_60Hz  = zeros(350,214);
fooof_offs = zeros(350,214);
fooof_offs_2_60Hz  = zeros(350,214);

for i = 1:350 
    fooof_exp(i,:) = all_subjects_exponent{i,1};
    fooof_exp_2_60Hz(i,:) = all_subjects_exponent_2_60Hz{i,1};
    fooof_offs(i,:) = all_subjects_offset{i,1};
    fooof_offs_2_60Hz(i,:) = all_subjects_offset_2_60Hz{i,1};
end

fooof_exp = fooof_exp(:,tmp_idx);
fooof_exp_2_60Hz = fooof_exp_2_60Hz(:,tmp_idx);
fooof_offs = fooof_offs(:,tmp_idx);
fooof_offs_2_60Hz = fooof_offs_2_60Hz(:,tmp_idx);

save([filepath '\DATA\foof_variables.mat'], "fooof_exp","fooof_offs_2_60Hz","fooof_offs","fooof_exp_2_60Hz");
writematrix(fooof_exp,[filepath '\DATA\fooof_exp.csv'])
writematrix(fooof_exp_2_60Hz,[filepath '\DATA\fooof_exp_2_60Hz.csv'])
writematrix(fooof_offs,[filepath '\DATA\fooof_offs.csv'])
writematrix(fooof_offs_2_60Hz,[filepath '\DATA\fooof_offs_2_60Hz.csv'])

%% FOR ADDITIONAL ANALYSIS: FIND Y-COORDINATES OF ALL PARCELS
% change the sorting of the parcel coordinates from fsaverage to schäfer200
% sort the rows of the averaged coordinates
centcoord_averaged = label_pos(tmp_idx,:);

% save coordinates
% save schäfer200 coordinates averaged in each parcel across vertices and
% sorted according to schäfer200
writematrix(centcoord_averaged,[pathin '\DATA\fsaverage_parcel_coorinates_averaged.csv'])
