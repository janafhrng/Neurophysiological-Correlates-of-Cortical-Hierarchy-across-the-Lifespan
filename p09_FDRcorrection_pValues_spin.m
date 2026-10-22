% FDR Correction HCTSA
clearvars;
close all
clc

filepath = pwd;

load([filepath '\DATA\lme_out_age_effect_hierarchy.mat']);

%% for the main effect of cortical hierarchy!
% FDR correction means --> p < 0.05/n-features
n_feat = 5961+7;
p_FDR = 0.05/n_feat;

% sort features for main effect of cortical hierarchy --> column 3
[~, tmp_idx] = sort(abs(allst(:,3)),'descend');

p_val_sort = pval(tmp_idx,3);

p_val_below = p_val_sort<p_FDR;

%% for the interaction effect of cortical hierarchy and age!
% FDR correction means --> p < 0.05/n-features
n_feat = 5961+7;
p_FDR = 0.05/n_feat;

% sort features for main effect of cortical hierarchy --> column 4
[~, tmp_idx] = sort(abs(allst(:,4)),'descend');

p_val_sort = pval(tmp_idx,4);

p_val_below = p_val_sort<p_FDR;




