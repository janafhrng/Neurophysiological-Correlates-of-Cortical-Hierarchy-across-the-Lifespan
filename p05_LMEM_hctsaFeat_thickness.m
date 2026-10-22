clearvars;
close all
clc

filepath = pwd;
 
%% data
% we decided to take participant number 306 out, becuase their cortical
% thickness value diviates too much from the others

%loading matrix with the data in the form of features x parcels
% matrix is already z-scored 
load([filepath '\DATA\mat_z.mat'])
mat_z(306,:,:) = [];

% loading the feature file with the complete list of features 
features = readtable([filepath '\DATA\features_z.csv']);

% loading cortical thickness
load([filepath '\DATA\thickness_smoothed_12_schaefer200_reordered.mat']);
thickness(306,:) = [];

demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));   
age_years(306) = [];

%% extract the statistics of the linear mixed effect model for each feature
for ifeature=1:size(mat_z,3) % loop over features
    if (mod(ifeature,200)==0)
        disp(ifeature);
    end

    subj=[];
    parcel=[];
    thick=[];
    feat=[];
    age=[];
    for k1=1:size(mat_z,1) %subjects
        subj=[subj ones(1,200)*k1];
        parcel=[parcel 1:200];
        thick=[thick thickness(k1,:)];
        feat=[feat squeeze(mat_z(k1,:,ifeature))];
        age = [age ones(1,200)*age_years(k1)];
    end
    varnames={'subj','parcel','feat','thick','age'};
    tbl=table(subj',parcel',feat',thick',age','VariableNames',varnames);
    model= 'feat ~  age*thick + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        %plotResiduals(lme,'fitted')
        allst(ifeature,:)=lme.Coefficients.tStat;
        pval(ifeature,:)=lme.Coefficients.pValue;
    end
end

% save lme statistics
save([filepath '\DATA\lme_hctsa_cortical_thickness_tstats.mat'],'allst','pval', '-v7.3');