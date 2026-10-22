% Test traditional time-series features and their relationship to cortica
% hierarchy

% time-series features: power in canonical freqbands;peak freq;1/f exponent
% LMEM with a: cortical hierarchy; b: cortical thickness

clearvars;
close all
clc

filepath = pwd;

% loading the the power of different freq bands over all parcels
load([filepath '\DATA\03_specData_freqbands_norm.mat'])

load([filepath '\DATA\hierarchy_parc.csv'])

demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));     

%% lme to find whether canonical power spectrum features can also descibe cortical hierachry?
freqbands{1} = zscore(delta_norm,[],2); freqbands{2} = zscore(theta_norm,[],2); freqbands{3} = zscore(alpha_norm,[],2); freqbands{4} = zscore(beta_norm,[],2); freqbands{5} = zscore(gamma_norm,[],2); 

% save for spin tests
for i = 1:5
    freqbands_zscored(i,:) = mean(freqbands{i},1);
end
writematrix(freqbands_zscored,[filepath '\DATA\03_specData_freqbands_norm_zscored.csv'])

% lmem (allst row1= intercept, row2=hierarchy, row3=age, row4=interation)
for ifreq = 1:5 % looping through the canonical frequency bands
    subj=[];
    parcel=[];
    freq=[];
    hierarchy=[];
    age=[];
    for k1=1:350 %subjects
        subj=[subj ones(1,200)*k1];
        freq=[freq freqbands{ifreq}(k1,:)];
        hierarchy=[hierarchy hierarchy_parc'];
        parcel = [parcel 1:200];
        age = [age ones(1,200)*age_years(k1)];
    end
    varnames={'subj','freq','hierarchy','parcel','age'};
    tbl=table(subj',freq',hierarchy',parcel',age','VariableNames',varnames);
    model= 'freq ~ age*hierarchy + (1|subj) + (1|parcel)';

    try
        lme = fitlme(tbl,model);
        allst_freqbands(ifreq,:)=lme.Coefficients.tStat;
        pval_freqbands(ifreq,:)=lme.Coefficients.pValue;
    end
end

%% SAME THING WITH CENTRE OF GRAVITY AND FOOOF EXPONENT
% loading center of gravity of the peak frequency
load([filepath '\DATA\peak_freq_cog.mat'])
peak_mean = mean(peak_freq_cog,1);

load([filepath '\DATA\foof_variables.mat'])

mat(:,:,1) = zscore(peak_freq_cog,[],2);
mat(:,:,2) = zscore(fooof_exp,[],2);
mat(:,:,3) = zscore(fooof_exp_2_60Hz,[],2);
mat(:,:,4) = zscore(fooof_offs,[],2);
mat(:,:,5) = zscore(fooof_offs_2_60Hz,[],2);

% extract the statistics of the linear mixed effect model for each feature
for ifeat = 1:5 % looping through the extra features 
    subj=[];
    parcel=[];
    feat=[];
    hierarchy=[];
    age=[];
    for k1=1:350 %subjects
        subj=[subj ones(1,200)*k1];
        feat=[feat mat(k1,:,ifeat)];
        hierarchy=[hierarchy hierarchy_parc'];
        parcel = [parcel 1:200];
        age = [age ones(1,200)*age_years(k1)];
    end
    varnames={'subj','feat','hierarchy','parcel','age'};
    tbl=table(subj',feat',hierarchy',parcel',age','VariableNames',varnames);
    model= 'feat ~ age*hierarchy + (1|subj) + (1|parcel)';

    try
        lme = fitlme(tbl,model);
        allst_extrafeat(ifeat,:)=lme.Coefficients.tStat;
        pval_extrafeat(ifeat,:)=lme.Coefficients.pValue;
    end
end





%% do the same with cortical thickness
% participant number 306 is an outlier related to corticla thickness, thus
% we do not use this data

for i = 1:5
    freqbands{i}(306,:) = [];
end

% loading cortical thickness
load([filepath '\DATA\thickness_smoothed_12_schaefer200_reordered.mat']);
thickness(306,:) = [];

age_years(306,:) = [];   

for ifreq = 1:5 % looping through the canonical frequency bands
    subj=[];
    parcel=[];
    freq=[];
    thick=[];
    age=[];
    for k1=1:349 %subjects
        subj=[subj ones(1,200)*k1];
        freq=[freq freqbands{ifreq}(k1,:)];
        parcel = [parcel 1:200];
        thick=[thick thickness(k1,:)];
        age = [age ones(1,200)*age_years(k1)];
    end
    varnames={'subj','freq','parcel','thick','age'};
    tbl=table(subj',freq',parcel',thick',age','VariableNames',varnames);
    model= 'freq ~ age*thick + (1|subj) + (1|parcel)';

    try
        lme = fitlme(tbl,model);
        allst_freqbands(ifreq,:)=lme.Coefficients.tStat;
        pval_freqbands(ifreq,:)=lme.Coefficients.pValue;
    end
end





%% SAME THING WITH CENTRE OF GRAVITY AND FOOOF EXPONENT AND CORTICAL THICKNESS
% loading center of gravity of the peak frequency

load([filepath '\DATA\foof_variables.mat'])
mat(:,:,1) = zscore(peak_freq_cog,[],2);
mat(:,:,2) = zscore(fooof_exp,[],2);
mat(:,:,3) = zscore(fooof_exp_2_60Hz,[],2);
mat(:,:,4) = zscore(fooof_offs,[],2);
mat(:,:,5) = zscore(fooof_offs_2_60Hz,[],2);
mat(306,:,:) = [];

% extract the statistics of the linear mixed effect model for each feature
for ifeat = 1:5 % looping through the extra features 
    subj=[];
    parcel=[];
    feat=[];
    thick=[];
    age=[];
    for k1=1:349 %subjects
        subj=[subj ones(1,200)*k1];
        feat=[feat mat(k1,:,ifeat)];
        thick=[thick thickness(k1,:)];
        parcel = [parcel 1:200];
        age = [age ones(1,200)*age_years(k1)];
    end
    varnames={'subj','feat','parcel','thick','age'};
    tbl=table(subj',feat',parcel',thick',age','VariableNames',varnames);
    model= 'feat ~ age*thick + (1|subj) + (1|parcel)';

    try
        lme = fitlme(tbl,model);
        allst_extrafeat(ifeat,:)=lme.Coefficients.tStat;
        pval_extrafeat(ifeat,:)=lme.Coefficients.pValue;
    end
end



