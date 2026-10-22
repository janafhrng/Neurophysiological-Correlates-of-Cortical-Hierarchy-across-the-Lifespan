% CONTROL ANALYSIS 
% control analysis for LMEM: spinning corticla hierarchy or corticla
% thickness to get a distribution of expectable t-values to test original t-value against

filepath = pwd;

% load the parcellated map of cortical hierarchy (sydnor et al., 2021)
load([filepath '\DATA\hierarchy_parc.csv'])

% load demogrpahics in order to catch the age in years
demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));      

% load Cortcial hierarchy permutations (50.000 according to alexander bloch)
permutations = readmatrix([filepath '\DATA\spins_hierarchy_permutations.csv']);

% load mat_z file
load('G:\JF\HCTSA feature gradients\CamCan\big DATA files\mat_z.mat')
features = readtable([filepath '\DATA\features_z.csv']);

%% AC at lag 133 ms
idx_feature = find(strcmp("AC_40",table2cell(features(:,"Var2"))));
mat_feat = squeeze(mat_z(:,:,idx_feature));

% extract the statistics of the linear mixed effect model for each feature
subj=[];
parcel=[];
age=[];
feat=[];
y=[];
for k1=1:350 %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    age = [age ones(1,200)*age_years(k1)];
    feat = [feat mat_feat(k1,:)]; 
end

for rot  = 1:10000
    perm=[];
    for s = 1:350
        perm =[perm permutations(:,rot)'];
    end

    varnames={'subj','parcel','feat','age','perm'};
    tbl=table(subj',parcel',feat',age',perm','VariableNames',varnames);
    model = 'feat ~ perm*age + (1|subj) + (1|parcel)';

    try
        lme=fitlme(tbl,model);

        %plotResiduals(lme,'fitted')
        lm(rot)=lme.Coefficients.tStat(3);
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(lm,"ascend");
lm_sorted(idxmax5) % 62.65
lm_sorted(idxmax1) % 80.2
lm_sorted(idxmin5) % -60.60
lm_sorted(idxmin1) % -74.67

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm,'NumBins',20)
xlabel('LMEM t-values for permutated AC40 map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];

hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)

xline(27.37,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100

box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_AC40_tvalues_2.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
p(1,1) =  (10000-max(find(lm_sorted<27.37)))/10000;


%% FREQUENCY BANDS!

% POWER SPECTRUM
% load power spectrum
load([filepath '\DATA\03_specData_freqbands_norm.mat']);

% z-score
delta_z = zscore(delta_norm,[],2);
theta_z = zscore(theta_norm,[],2);
alpha_z = zscore(alpha_norm,[],2);
beta_z = zscore(beta_norm,[],2);
gamma_z = zscore(gamma_norm,[],2);

% additional variables
load([filepath '\DATA\peak_freq_cog.mat'])
load([filepath '\DATA\foof_variables.mat'])
cog = zscore(peak_freq_cog,[],2);
fooof = zscore(fooof_exp_2_60Hz,[],2);


%% DELTA

% extract the statistics of the linear mixed effect model for each feature
subj=[];
parcel=[];
age=[];
feat=[];
for k1=1:350 %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    age = [age ones(1,200)*age_years(k1)];
    feat = [feat squeeze(delta_z(k1,:))];
end

for rot  = 1261:10000
    if ~mod(rot,500)
        rot
    end

    perm=[];
    for s = 1:350
        perm =[perm permutations(:,rot)'];
    end

    varnames={'subj','parcel','feat','age','perm'};
    tbl=table(subj',parcel',feat',age',perm','VariableNames',varnames);
    model = 'feat ~ perm*age + (1|subj) + (1|parcel)';

    try
        lme=fitlme(tbl,model);

        %plotResiduals(lme,'fitted')
        lm(rot)=lme.Coefficients.tStat(3);
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(lm,"ascend");
lm_sorted(idxmax5) 
lm_sorted(idxmax1)
lm_sorted(idxmin5) 
lm_sorted(idxmin1) 

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm,'NumBins',20)
xlabel('LMEM t-values for permutated delta map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];

hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)

xline(16.88,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100

box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_delta_tvalues_2.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
p(1,2) =  (10000-max(find(lm_sorted<16.88)))/10000;

%% THETA
lm =[];

% extract the statistics of the linear mixed effect model for each feature
subj=[];
parcel=[];
age=[];
feat=[];
for k1=1:350 %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    age = [age ones(1,200)*age_years(k1)];
    feat = [feat squeeze(theta_z(k1,:))];
end

for rot  = 1:10000
    if ~mod(rot,500)
        rot
    end


% which p-value does the t-vual correspond to?
p(1,2) =  (10000-max(find(lm_sorted<16.88)))/10000;
    perm=[];
    for s = 1:350
        perm =[perm permutations(:,rot)'];
    end

    varnames={'subj','parcel','feat','age','perm'};
    tbl=table(subj',parcel',feat',age',perm','VariableNames',varnames);
    model = 'feat ~ perm*age + (1|subj) + (1|parcel)';

    try
        lme=fitlme(tbl,model);

        %plotResiduals(lme,'fitted')
        lm(rot)=lme.Coefficients.tStat(3);
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(lm,"ascend");
lm_sorted(idxmax5) 
lm_sorted(idxmax1) 
lm_sorted(idxmin5) 
lm_sorted(idxmin1)

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm,'NumBins',20)
xlabel('LMEM t-values for permutated theta map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];

hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)

xline(19.08,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100

box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_theta_tvalues_2.png'],'Resolution',500)


% which p-value does the t-vual correspond to?
p(1,3) =  (10000-max(find(lm_sorted<19.08)))/10000;


%% ALPHA

% extract the statistics of the linear mixed effect model for each feature
subj=[];
parcel=[];
age=[];
feat=[];
for k1=1:350 %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    age = [age ones(1,200)*age_years(k1)];
    feat = [feat squeeze(alpha_z(k1,:))];
end

for rot  = 1:10000
    if ~mod(rot,500)
        rot
    end

    perm=[];
    for s = 1:350
        perm =[perm permutations(:,rot)'];
    end

    varnames={'subj','parcel','feat','age','perm'};
    tbl=table(subj',parcel',feat',age',perm','VariableNames',varnames);
    model = 'feat ~ perm*age + (1|subj) + (1|parcel)';

    try
        lme=fitlme(tbl,model);

        %plotResiduals(lme,'fitted')
        lm(rot)=lme.Coefficients.tStat(3);
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(lm,"ascend");
lm_sorted(idxmax5) 
lm_sorted(idxmax1) 
lm_sorted(idxmin5) 
lm_sorted(idxmin1) 

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm,'NumBins',20)
xlabel('LMEM t-values for permutated alpha map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];

hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)

xline(-12.31,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100

box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_alpha_tvalues_2.png'],'Resolution',500)


% which p-value does the t-vual correspond to?
p(1,4) =  (10000-max(find(lm_sorted<-12.31)))/10000;



%% BETA 

% extract the statistics of the linear mixed effect model for each feature
subj=[];
parcel=[];
age=[];
feat=[];
for k1=1:350 %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    age = [age ones(1,200)*age_years(k1)];
    feat = [feat squeeze(beta_z(k1,:))];
end

for rot  = 1:10000
    if ~mod(rot,500)
        rot
    end

    perm=[];
    for s = 1:350
        perm =[perm permutations(:,rot)'];
    end

    varnames={'subj','parcel','feat','age','perm'};
    tbl=table(subj',parcel',feat',age',perm','VariableNames',varnames);
    model = 'feat ~ perm*age + (1|subj) + (1|parcel)';

    try
        lme=fitlme(tbl,model);

        %plotResiduals(lme,'fitted')
        lm(rot)=lme.Coefficients.tStat(3);
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(lm,"ascend");
lm_sorted(idxmax5) 
lm_sorted(idxmax1) 
lm_sorted(idxmin5) 
lm_sorted(idxmin1) 

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm,'NumBins',20)  
xlabel('LMEM t-values for permutated beta map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];

hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)

xline(-0.83,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100

box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_beta_tvalues_2.png'],'Resolution',500)



% which p-value does the t-vual correspond to?
p(1,5) =  (10000-max(find(lm_sorted<-0.83)))/10000;


%% GAMMA 

% extract the statistics of the linear mixed effect model for each feature
subj=[];
parcel=[];
age=[];
feat=[];
for k1=1:350 %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    age = [age ones(1,200)*age_years(k1)];
    feat = [feat squeeze(gamma_z(k1,:))];
end

for rot  = 1:10000
    if ~mod(rot,500)
        rot
    end

    perm=[];
    for s = 1:350
        perm =[perm permutations(:,rot)'];
    end

    varnames={'subj','parcel','feat','age','perm'};
    tbl=table(subj',parcel',feat',age',perm','VariableNames',varnames);
    model = 'feat ~ perm*age + (1|subj) + (1|parcel)';

    try
        lme=fitlme(tbl,model);

        %plotResiduals(lme,'fitted')
        lm(rot)=lme.Coefficients.tStat(3);
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(lm,"ascend");
lm_sorted(idxmax5) 
lm_sorted(idxmax1) 
lm_sorted(idxmin5) 
lm_sorted(idxmin1)

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm,'NumBins',20)  
xlabel('LMEM t-values for permutated gamma map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];

hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)

xline(12.7,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100

box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_gamma_tvalues_2.png'],'Resolution',500)



% which p-value does the t-vual correspond to?
p(1,6) =  (10000-max(find(lm_sorted<12.7)))/10000;



%% PEAK FREQUENCY

% extract the statistics of the linear mixed effect model for each feature
subj=[];
parcel=[];
age=[];
feat=[];
for k1=1:350 %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    age = [age ones(1,200)*age_years(k1)];
    feat = [feat squeeze(cog(k1,:))];
end

for rot  = 1:10000
    if ~mod(rot,500)
        rot
    end

    perm=[];
    for s = 1:350
        perm =[perm permutations(:,rot)'];
    end

    varnames={'subj','parcel','feat','age','perm'};
    tbl=table(subj',parcel',feat',age',perm','VariableNames',varnames);
    model = 'feat ~ perm*age + (1|subj) + (1|parcel)';

    try
        lme=fitlme(tbl,model);

        %plotResiduals(lme,'fitted')
        lm(rot)=lme.Coefficients.tStat(3);
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(lm,"ascend");
lm_sorted(idxmax5) 
lm_sorted(idxmax1) 
lm_sorted(idxmin5) 
lm_sorted(idxmin1) 

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm,'NumBins',20)
xlabel('LMEM t-values for permutated cog map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];

hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)

xline(-10.34,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100

box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_cog_tvalues_2.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
p(1,7) =  (10000-max(find(lm_sorted<-10.34)))/10000;



% %% FOOOF
% 
% % extract the statistics of the linear mixed effect model for each feature
% subj=[];
% parcel=[];
% age=[];
% feat=[];
% for k1=1:350 %subjects
%     subj=[subj ones(1,200)*k1];
%     parcel=[parcel 1:200];
%     age = [age ones(1,200)*age_years(k1)];
%     feat = [feat squeeze(fooof(k1,:))];
% end
% 
% for rot  = 1:10000
%     perm=[];
%     for s = 1:350
%         perm =[perm permutations(:,rot)'];
%     end
% 
%     varnames={'subj','parcel','feat','age','perm'};
%     tbl=table(subj',parcel',feat',age',perm','VariableNames',varnames);
%     model = 'perm ~ feat*age + (1|subj) + (1|parcel)';
% 
%     try
%         lme=fitlme(tbl,model);
% 
%         %plotResiduals(lme,'fitted')
%         lm(rot)=lme.Coefficients.tStat(2);
%     end
% end
% 
% % find indices for 2.5% and 97.5% where the t-valzes should be significant
% idxmin5 = 0.025*10000;
% idxmin1 = 0.005*10000;
% idxmax5 = 0.975*10000;
% idxmax1 = 0.995*10000;
% 
% lm_sorted = sort(lm,"ascend");
% lm_sorted(idxmax5) % 62.65
% lm_sorted(idxmax1) % 80.2
% lm_sorted(idxmin5) % -60.60
% lm_sorted(idxmin1) % -74.67
% 
% figure('Color','w', 'Position', [0 0 600 500])
% h = histogram(lm,'NumBins',20)
% xlabel('LMEM t-values for permutated fooof map')
% ylabel('count')
% h.FaceColor = [0 0.5 0.5];
% 
% hold on 
% xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
% 
% xline(-33.2,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
% 
% box off
% exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_fooof_tvalues_2.png'],'Resolution',500)
% 
% % which p-value does the t-vual correspond to?
% p(1,8) =  (10000-max(find(lm_sorted<-33.2)))/10000;
% 









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%% CORTICAL THICKNESS
% do thesame with cortical thickness
spin_index = readmatrix([filepath, '\DATA\spins_index_permutations.csv']);

% loading cortical thickness
load([filepath '\DATA\thickness_smoothed_12_schaefer200_reordered.mat']);
thickness(306,:) = [];

demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));     
age_years(306,:) = [];   

% loading the the power of different freq bands over all parcels
load([filepath '\DATA\03_specData_freqbands_norm.mat'])
% z-score
delta_z = zscore(delta_norm,[],2); delta_z(306,:) = [];
theta_z = zscore(theta_norm,[],2); theta_z(306,:) = [];
alpha_z = zscore(alpha_norm,[],2); alpha_z(306,:) = [];
beta_z = zscore(beta_norm,[],2); beta_z(306,:) = [];
gamma_z = zscore(gamma_norm,[],2); gamma_z(306,:) = [];

% additional variables
load([filepath '\DATA\peak_freq_cog.mat'])
load([filepath '\DATA\foof_variables.mat'])
cog = zscore(peak_freq_cog,[],2); cog(306,:) = [];
fooof = zscore(fooof_exp_2_60Hz,[],2); fooof(306,:) = [];

% AC at lag 133 ms
idx_feature = find(strcmp("AC_40",table2cell(features(:,"Var2"))));
mat_feat = squeeze(mat_z(:,:,idx_feature));
mat_feat(306,:) = [];


%% AC40
% extract the statistics of the linear mixed effect model for each feature
subj=[];
parcel=[];
age=[];
feat=[];
for k1=1:349 %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    age = [age ones(1,200)*age_years(k1)];
    feat = [feat mat_feat(k1,:)]; 
end

for rot  = 1:10000
    perm=[];
    for s = 1:349
        perm =[perm thickness(s,spin_index(:,rot))];
    end

    varnames={'subj','parcel','feat','age','perm'};
    tbl=table(subj',parcel',feat',age',perm','VariableNames',varnames);
    model = 'perm ~ feat*age + (1|subj) + (1|parcel)';

    try
        lme=fitlme(tbl,model);

        %plotResiduals(lme,'fitted')
        lm(rot)=lme.Coefficients.tStat(2);
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(lm,"ascend");
lm_sorted(idxmax5) % 62.65
lm_sorted(idxmax1) % 80.2
lm_sorted(idxmin5) % -60.60
lm_sorted(idxmin1) % -74.67

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm,'NumBins',20)
xlabel('LMEM t-values for permutated AC40 map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];

hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)

xline(19.43,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100

box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_AC40_tvalues_thick_2.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
p(2,1) =  (10000-max(find(lm_sorted<19.43)))/10000;


%% DELTA

% extract the statistics of the linear mixed effect model for each feature
subj=[];
parcel=[];
age=[];
feat=[];
for k1=1:349 %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    age = [age ones(1,200)*age_years(k1)];
    feat = [feat squeeze(delta_z(k1,:))];
end

for rot  = 1:10000
    if ~mod(rot,500)
        rot
    end

    perm=[];
    for s = 1:349
        perm =[perm thickness(s,spin_index(:,rot))];
    end

    varnames={'subj','parcel','feat','age','perm'};
    tbl=table(subj',parcel',feat',age',perm','VariableNames',varnames);
    model = 'feat ~ perm*age + (1|subj) + (1|parcel)';

    try
        lme=fitlme(tbl,model);

        %plotResiduals(lme,'fitted')
        lm(rot)=lme.Coefficients.tStat(3);
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(lm,"ascend");
lm_sorted(idxmax5) 
lm_sorted(idxmax1)
lm_sorted(idxmin5) 
lm_sorted(idxmin1) 

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm,'NumBins',20)
xlabel('LMEM t-values for permutated delta map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];

hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)

xline(23.86,'r--','LineWidth',2,'Alpha',0.7) % delta

box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_delta_tvalues_thick_2.png'],'Resolution',500)


% which p-value does the t-vual correspond to?
p(2,2) =  (10000-max(find(lm_sorted<23.86)))/10000;


%% THETA

% extract the statistics of the linear mixed effect model for each feature
% subj=[];
% parcel=[];
% age=[];
feat=[];
for k1=1:349 %subjects
%     subj=[subj ones(1,200)*k1];
%     parcel=[parcel 1:200];
%     age = [age ones(1,200)*age_years(k1)];
    feat = [feat squeeze(theta_z(k1,:))];
end

for rot  = 1:10000
    if ~mod(rot,500)
        rot
    end

    perm=[];
    for s = 1:349
        perm =[perm thickness(s,spin_index(:,rot))];
    end

    varnames={'subj','parcel','feat','age','perm'};
    tbl=table(subj',parcel',feat',age',perm','VariableNames',varnames);
    model = 'feat ~ perm*age + (1|subj) + (1|parcel)';

    try
        lme=fitlme(tbl,model);

        %plotResiduals(lme,'fitted')
        lm(rot)=lme.Coefficients.tStat(3);
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(lm,"ascend");
lm_sorted(idxmax5) 
lm_sorted(idxmax1) 
lm_sorted(idxmin5) 
lm_sorted(idxmin1) 

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm,'NumBins',20)
xlabel('LMEM t-values for permutated theta map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];

hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)

xline(16.98,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100

box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_theta_tvalues_thick_2.png'],'Resolution',500)


% which p-value does the t-vual correspond to?
p(2,3) =  (10000-max(find(lm_sorted<16.98)))/10000;


%% ALPHA

% extract the statistics of the linear mixed effect model for each feature
% subj=[];
% parcel=[];
% age=[];
feat=[];
for k1=1:349 %subjects
%     subj=[subj ones(1,200)*k1];
%     parcel=[parcel 1:200];
%     age = [age ones(1,200)*age_years(k1)];
    feat = [feat squeeze(alpha_z(k1,:))];
end

for rot  = 1:10000
    if ~mod(rot,500)
        rot
    end

    perm=[];
    for s = 1:349
        perm =[perm thickness(s,spin_index(:,rot))];
    end

    varnames={'subj','parcel','feat','age','perm'};
    tbl=table(subj',parcel',feat',age',perm','VariableNames',varnames);
    model = 'feat ~ perm*age + (1|subj) + (1|parcel)';

    try
        lme=fitlme(tbl,model);

        %plotResiduals(lme,'fitted')
        lm(rot)=lme.Coefficients.tStat(3);
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(lm,"ascend");
lm_sorted(idxmax5) % 62.65
lm_sorted(idxmax1) % 80.2
lm_sorted(idxmin5) % -60.60
lm_sorted(idxmin1) % -74.67

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm,'NumBins',20)
xlabel('LMEM t-values for permutated alpha map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];

hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)

xline(-29.47,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100

box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_alpha_tvalues_thick_2.png'],'Resolution',500)


% which p-value does the t-vual correspond to?
p(2,4) =  (10000-max(find(lm_sorted<-29.47)))/10000;



%% BETA

% extract the statistics of the linear mixed effect model for each feature
% subj=[];
% parcel=[];
% age=[];
feat=[];
for k1=1:349 %subjects
%     subj=[subj ones(1,200)*k1];
%     parcel=[parcel 1:200];
%     age = [age ones(1,200)*age_years(k1)];
    feat = [feat squeeze(beta_z(k1,:))];
end

for rot  = 1:10000
    if ~mod(rot,500)
        rot
    end

    perm=[];
    for s = 1:349
        perm =[perm thickness(s,spin_index(:,rot))];
    end

    varnames={'subj','parcel','feat','age','perm'};
    tbl=table(subj',parcel',feat',age',perm','VariableNames',varnames);
    model = 'feat ~ perm*age + (1|subj) + (1|parcel)';

    try
        lme=fitlme(tbl,model);

        %plotResiduals(lme,'fitted')
        lm(rot)=lme.Coefficients.tStat(3);
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(lm,"ascend");
lm_sorted(idxmax5) % 62.65
lm_sorted(idxmax1) % 80.2
lm_sorted(idxmin5) % -60.60
lm_sorted(idxmin1) % -74.67

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm,'NumBins',20)
xlabel('LMEM t-values for permutated beta map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];

hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)

xline(-0.25,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100

box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_beta_tvalues_thick_2.png'],'Resolution',500)


% which p-value does the t-vual correspond to?
p(2,5) =  (10000-max(find(lm_sorted<-0.25)))/10000;



%% GAMMA

% extract the statistics of the linear mixed effect model for each feature
% subj=[];
% parcel=[];
% age=[];
feat=[];
for k1=1:349 %subjects
%     subj=[subj ones(1,200)*k1];
%     parcel=[parcel 1:200];
%     age = [age ones(1,200)*age_years(k1)];
    feat = [feat squeeze(gamma_z(k1,:))];
end

for rot  = 1:10000
    if ~mod(rot,500)
        rot
    end

    perm=[];
    for s = 1:349
        perm =[perm thickness(s,spin_index(:,rot))];
    end

    varnames={'subj','parcel','feat','age','perm'};
    tbl=table(subj',parcel',feat',age',perm','VariableNames',varnames);
    model = 'feat ~ perm*age + (1|subj) + (1|parcel)';

    try
        lme=fitlme(tbl,model);

        %plotResiduals(lme,'fitted')
        lm(rot)=lme.Coefficients.tStat(3);
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(lm,"ascend");
lm_sorted(idxmax5) % 62.65
lm_sorted(idxmax1) % 80.2
lm_sorted(idxmin5) % -60.60
lm_sorted(idxmin1) % -74.67

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm,'NumBins',20)
xlabel('LMEM t-values for permutated gamma map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];

hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)

xline(18,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100

box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_gamma_tvalues_thick_2.png'],'Resolution',500)



% which p-value does the t-vual correspond to?
p(2,6) =  (10000-max(find(lm_sorted<18)))/10000;




%% COG

% extract the statistics of the linear mixed effect model for each feature
% subj=[];
% parcel=[];
% age=[];
feat=[];
for k1=1:349 %subjects
%     subj=[subj ones(1,200)*k1];
%     parcel=[parcel 1:200];
%     age = [age ones(1,200)*age_years(k1)];
    feat = [feat squeeze(cog(k1,:))];
end

for rot  = 1:10000
    if ~mod(rot,500)
        rot
    end

    perm=[];
    for s = 1:349
        perm =[perm thickness(s,spin_index(:,rot))];
    end

    varnames={'subj','parcel','feat','age','perm'};
    tbl=table(subj',parcel',feat',age',perm','VariableNames',varnames);
    model = 'feat ~ perm*age + (1|subj) + (1|parcel)';

    try
        lme=fitlme(tbl,model);

        %plotResiduals(lme,'fitted')
        lm(rot)=lme.Coefficients.tStat(3);
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(lm,"ascend");
lm_sorted(idxmax5) % 62.65
lm_sorted(idxmax1) % 80.2
lm_sorted(idxmin5) % -60.60
lm_sorted(idxmin1) % -74.67

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm,'NumBins',20)
xlabel('LMEM t-values for permutated cog map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];

hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)

xline(-7.09,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100

box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_cog_tvalues_thick_2.png'],'Resolution',500)



% which p-value does the t-vual correspond to?
p(2,7) =  (10000-max(find(lm_sorted<-7.09)))/10000;



% %% FOOOF
% 
% % extract the statistics of the linear mixed effect model for each feature
% % subj=[];
% % parcel=[];
% % age=[];
% feat=[];
% for k1=1:349 %subjects
%     %     subj=[subj ones(1,200)*k1];
%     %     parcel=[parcel 1:200];
%     %     age = [age ones(1,200)*age_years(k1)];
%     feat = [feat squeeze(fooof(k1,:))];
% end
% 
% for rot  = 1:10000
%     perm=[];
%     for s = 1:349
%         perm =[perm thickness(s,spin_index(:,rot))];
%     end
% 
%     varnames={'subj','parcel','feat','age','perm'};
%     tbl=table(subj',parcel',feat',age',perm','VariableNames',varnames);
%     model = 'perm ~ feat*age + (1|subj) + (1|parcel)';
% 
%     try
%         lme=fitlme(tbl,model);
% 
%         %plotResiduals(lme,'fitted')
%         lm(rot)=lme.Coefficients.tStat(2);
%     end
% end
% 
% % find indices for 2.5% and 97.5% where the t-valzes should be significant
% idxmin5 = 0.025*10000;
% idxmin1 = 0.005*10000;
% idxmax5 = 0.975*10000;
% idxmax1 = 0.995*10000;
% 
% lm_sorted = sort(lm,"ascend");
% lm_sorted(idxmax5) % 62.65
% lm_sorted(idxmax1) % 80.2
% lm_sorted(idxmin5) % -60.60
% lm_sorted(idxmin1) % -74.67
% 
% figure('Color','w', 'Position', [0 0 600 500])
% h = histogram(lm,'NumBins',20)
% xlabel('LMEM t-values for permutated fooof map')
% ylabel('count')
% h.FaceColor = [0 0.5 0.5];
% 
% hold on
% xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
% 
% xline(-16.2,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
% 
% box off
% exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_fooof_tvalues_thick_2.png'],'Resolution',500)
% 
% 
% 
% % which p-value does the t-vual correspond to?
% p(2,8) =  (10000-max(find(lm_sorted<-16.2)))/10000;
% 
% 
% 






















% %% what kind of t-values can be expected?
%
% filepath = pwd;
%
% %% STATAVL100
% % use the permutations from python where we create 50.000 permutations
% % with alexander bloch
% permutations = readmatrix([filepath '\DATA\spins_StatAvl100_permutations.csv']);
%
% % load demogrpahics in order to catch the age in years
% demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
% age_years = table2array(demographics(:,2));
%
% % load the parcellated map of cortical hierarchy (sydnor et al., 2021)
% load([filepath '\DATA\hierarchy_parc.csv'])
%
%
% %% null distribution for randomly expactable t-values
%
%
% % extract the statistics of the linear mixed effect model for each feature
% subj=[];
% parcel=[];
% age=[];
% hierarchy=[];
% y=[];
% for k1=1:350 %subjects
%     subj=[subj ones(1,200)*k1];
%     parcel=[parcel 1:200];
%     age = [age ones(1,200)*age_years(k1)];
%     hierarchy=[hierarchy hierarchy_parc'];
% end
%
%
% for rot  = 1:50000
%     perm=[];
%     for s = 1:350
%         perm =[perm permutations(:,rot)'];
%     end
%
%     varnames={'subj','parcel','perm','age','hierarchy'};
%     tbl=table(subj',parcel',perm',age',hierarchy','VariableNames',varnames);
%     model = 'perm ~ hierarchy*age + (1|subj)';
%
%     try
%         lme=fitlme(tbl,model);
%
%         %plotResiduals(lme,'fitted')
%         lm(rot)=lme.Coefficients.tStat(2);
%     end
% end
%
% % find indices for 2.5% and 97.5% where the t-valzes should be significant
% idxmin5 = 0.025*50000;
% idxmin1 = 0.005*50000;
% idxmax5 = 0.975*50000;
% idxmax1 = 0.995*50000;
%
% lm_sorted = sort(lm,"ascend");
% lm_sorted(idxmax5) % 62.65
% lm_sorted(idxmax1) % 80.2
% lm_sorted(idxmin5) % -60.60
% lm_sorted(idxmin1) % -74.67
%
% figure('Color','w', 'Position', [0 0 600 500])
% h = histogram(lm,'NumBins',20)
% xlabel('LMEM t-values for permutated StatAvl100 map')
% ylabel('count')
% h.FaceColor = [0 0.5 0.5];
%
% hold on
% xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
%
% xline(80.11,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
%
% box off
%
% exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_statavl_tvalues.png'],'Resolution',500)
%
%
%
% %% DELTA
% % use the permutations from python where we create 50.000 permutations
% % with alexander bloch
% permutations = readmatrix([filepath '\DATA\spins_delta_permutations.csv']);
%
% % extract the statistics of the linear mixed effect model for each feature
% subj=[];
% parcel=[];
% age=[];
% hierarchy=[];
% y=[];
% for k1=1:350 %subjects
%     subj=[subj ones(1,200)*k1];
%     parcel=[parcel 1:200];
%     age = [age ones(1,200)*age_years(k1)];
%     hierarchy=[hierarchy hierarchy_parc'];
% end
%
%
% for rot  = 1:50000
%     perm=[];
%     for s = 1:350
%         perm =[perm permutations(:,rot)'];
%     end
% 
%     varnames={'subj','parcel','perm','age','hierarchy'};
%     tbl=table(subj',parcel',perm',age',hierarchy','VariableNames',varnames);
%     model = 'perm ~ hierarchy*age + (1|subj)';
% 
%     try
%         lme=fitlme(tbl,model);
% 
%         %plotResiduals(lme,'fitted')
%         lm2(rot)=lme.Coefficients.tStat(3);
%     end
% end
% 
% % find indices for 2.5% and 97.5% where the t-valzes should be significant
% idxmin5 = 0.025*50000;
% idxmin1 = 0.005*50000;
% idxmax5 = 0.975*50000;
% idxmax1 = 0.995*50000;
% 
% lm_sorted = sort(lm2,"ascend");
% lm_sorted(idxmax5) % 62.65
% lm_sorted(idxmax1) % 80.2
% lm_sorted(idxmin5) % -60.60
% lm_sorted(idxmin1) % -74.67
% 
% figure('Color','w', 'Position', [0 0 600 500])
% h = histogram(lm2,'NumBins',20)
% xlabel('LMEM t-values for permutated delta power map')
% ylabel('count')
% h.FaceColor = [0 0.5 0.5];
% 
% hold on 
% xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
% 
% xline(78.1,'r--','LineWidth',2,'Alpha',0.7) % Delta
% 
% box off
% 
% exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_delta_tvalues.png'],'Resolution',500)
% 
% 
% %% THETA
% % use the permutations from python where we create 50.000 permutations
% % with alexander bloch 
% permutations = readmatrix([filepath '\DATA\spins_theta_permutations.csv']);
% 
% % extract the statistics of the linear mixed effect model for each feature
% subj=[];
% parcel=[];
% age=[];
% hierarchy=[];
% y=[];
% for k1=1:350 %subjects
%     subj=[subj ones(1,200)*k1];
%     parcel=[parcel 1:200];
%     age = [age ones(1,200)*age_years(k1)];
%     hierarchy=[hierarchy hierarchy_parc'];
% end
% 
% 
% for rot  = 1:50000
%     perm=[];
%     for s = 1:350
%         perm =[perm permutations(:,rot)'];
%     end
% 
%     varnames={'subj','parcel','perm','age','hierarchy'};
%     tbl=table(subj',parcel',perm',age',hierarchy','VariableNames',varnames);
%     model = 'perm ~ hierarchy*age + (1|subj)';
% 
%     try
%         lme=fitlme(tbl,model);
% 
%         %plotResiduals(lme,'fitted')
%         lm3(rot)=lme.Coefficients.tStat(3);
%     end
% end
% 
% % find indices for 2.5% and 97.5% where the t-valzes should be significant
% idxmin5 = 0.025*50000;
% idxmin1 = 0.005*50000;
% idxmax5 = 0.975*50000;
% idxmax1 = 0.995*50000;
% 
% lm_sorted = sort(lm3,"ascend");
% lm_sorted(idxmax5)
% lm_sorted(idxmax1) 
% lm_sorted(idxmin5)
% lm_sorted(idxmin1)
% 
% figure('Color','w', 'Position', [0 0 600 500])
% h = histogram(lm3,'NumBins',20)
% xlabel('LMEM t-values for permutated theta power map')
% ylabel('count')
% h.FaceColor = [0 0.5 0.5];
% 
% hold on 
% xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
% 
% xline(50.4,'r--','LineWidth',2,'Alpha',0.7) % t-value for theta power
% 
% box off
% 
% exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_theta_tvalues.png'],'Resolution',500)
% 
% 
% 
% %% ALPHA
% % use the permutations from python where we create 50.000 permutations
% % with alexander bloch 
% permutations = readmatrix([filepath '\DATA\spins_alpha_permutations.csv']);
% 
% % extract the statistics of the linear mixed effect model for each feature
% subj=[];
% parcel=[];
% age=[];
% hierarchy=[];
% y=[];
% for k1=1:350 %subjects
%     subj=[subj ones(1,200)*k1];
%     parcel=[parcel 1:200];
%     age = [age ones(1,200)*age_years(k1)];
%     hierarchy=[hierarchy hierarchy_parc'];
% end
% 
% 
% for rot  = 1:50000
%     perm=[];
%     for s = 1:350
%         perm =[perm permutations(:,rot)'];
%     end
% 
%     varnames={'subj','parcel','perm','age','hierarchy'};
%     tbl=table(subj',parcel',perm',age',hierarchy','VariableNames',varnames);
%     model = 'perm ~ hierarchy*age + (1|subj)';
% 
%     try
%         lme=fitlme(tbl,model);
% 
%         %plotResiduals(lme,'fitted')
%         lm4(rot)=lme.Coefficients.tStat(3);
%     end
% end
% 
% % find indices for 2.5% and 97.5% where the t-valzes should be significant
% idxmin5 = 0.025*50000;
% idxmin1 = 0.005*50000;
% idxmax5 = 0.975*50000;
% idxmax1 = 0.995*50000;
% 
% lm_sorted = sort(lm4,"ascend");
% lm_sorted(idxmax5)
% lm_sorted(idxmax1) 
% lm_sorted(idxmin5)
% lm_sorted(idxmin1)
% 
% figure('Color','w', 'Position', [0 0 600 500])
% h = histogram(lm4,'NumBins',20)
% xlabel('LMEM t-values for permutated alpha power map')
% ylabel('count')
% h.FaceColor = [0 0.5 0.5];
% 
% hold on 
% xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
% 
% xline(65.6,'r--','LineWidth',2,'Alpha',0.7) % t-value for alpha power
% 
% box off
% 
% exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_alpha_tvalues.png'],'Resolution',500)
% 
% 
% 
% %% GAMMA
% % use the permutations from python where we create 50.000 permutations
% % with alexander bloch 
% permutations = readmatrix([filepath '\DATA\spins_gamma_permutations.csv']);
% 
% % extract the statistics of the linear mixed effect model for each feature
% subj=[];
% parcel=[];
% age=[];
% hierarchy=[];
% y=[];
% for k1=1:350 %subjects
%     subj=[subj ones(1,200)*k1];
%     parcel=[parcel 1:200];
%     age = [age ones(1,200)*age_years(k1)];
%     hierarchy=[hierarchy hierarchy_parc'];
% end
% 
% 
% for rot  = 1:50000
%     perm=[];
%     for s = 1:350
%         perm =[perm permutations(:,rot)'];
%     end
% 
%     varnames={'subj','parcel','perm','age','hierarchy'};
%     tbl=table(subj',parcel',perm',age',hierarchy','VariableNames',varnames);
%     model = 'perm ~ hierarchy*age + (1|subj)';
% 
%     try
%         lme=fitlme(tbl,model);
% 
%         %plotResiduals(lme,'fitted')
%         lm5(rot)=lme.Coefficients.tStat(3);
%     end
% end
% 
% % find indices for 2.5% and 97.5% where the t-valzes should be significant
% idxmin5 = 0.025*50000;
% idxmin1 = 0.005*50000;
% idxmax5 = 0.975*50000;
% idxmax1 = 0.995*50000;
% 
% lm_sorted = sort(lm5,"ascend");
% lm_sorted(idxmax5) 
% lm_sorted(idxmax1) 
% lm_sorted(idxmin5)
% lm_sorted(idxmin1)
% 
% figure('Color','w', 'Position', [0 0 600 500])
% h = histogram(lm5,'NumBins',20)
% xlabel('LMEM t-values for permutated gamma power map')
% ylabel('count')
% h.FaceColor = [0 0.5 0.5];
% 
% hold on 
% xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
% 
% xline(64.9,'r--','LineWidth',2,'Alpha',0.7) % t-value for gamma power
% 
% box off
% 
% exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_gamma_tvalues.png'],'Resolution',500)










% 
% 
% %%%%%%%%%%%%%%%%%%%% DO THE SAME CONTROL ONLY WITH CORTICAL THICKNESS
% 
% 
% 
% 
% 
% 
% % loading cortical thickness
% load([filepath '\DATA\thickness_smoothed_12_schaefer200_reordered.mat']);
% thickness(306,:) = [];
% thickness = mean(thickness,1);
% 
% 
% % STATAVL100 % null distribution for randomly expactable t-values
% 
% % use the permutations from python where we create 50.000 permutations
% % with alexander bloch 
% permutations = readmatrix([filepath '\DATA\spins_StatAvl100_permutations.csv']);
% 
% 
% % extract the statistics of the linear mixed effect model for each feature
% subj=[];
% parcel=[];
% age=[];
% thick=[];
% y=[];
% for k1=1:350 %subjects
%     subj=[subj ones(1,200)*k1];
%     parcel=[parcel 1:200];
%     age = [age ones(1,200)*age_years(k1)];
%     thick=[thick thickness];
% end
% 
% 
% for rot  = 1:50000
%     perm=[];
%     for s = 1:350
%         perm =[perm permutations(:,rot)'];
%     end
% 
%     varnames={'subj','parcel','perm','age','thick'};
%     tbl=table(subj',parcel',perm',age',thick','VariableNames',varnames);
%     model = 'perm ~ thick*age + (1|subj)';
% 
%     try
%         lme=fitlme(tbl,model);
% 
%         %plotResiduals(lme,'fitted')
%         lm(rot)=lme.Coefficients.tStat(3);
%     end
% end
% 
% % find indices for 2.5% and 97.5% where the t-valzes should be significant
% idxmin5 = 0.025*50000;
% idxmin1 = 0.005*50000;
% idxmax5 = 0.975*50000;
% idxmax1 = 0.995*50000;
% 
% lm_sorted = sort(lm,"ascend");
% lm_sorted(idxmax5) % 62.65
% lm_sorted(idxmax1) % 80.2
% lm_sorted(idxmin5) % -60.60
% lm_sorted(idxmin1) % -74.67
% 
% figure('Color','w', 'Position', [0 0 600 500])
% h = histogram(lm,'NumBins',20)
% xlabel('LMEM t-values for permutated StatAvl100 map')
% ylabel('count')
% h.FaceColor = [0 0.5 0.5];
% 
% hold on 
% xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
% 
% xline(56.6,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
% 
% box off
% 
% exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_statavl_tvalues_thick.png'],'Resolution',500)
% 
% 
% 
% %% DELTA
% % use the permutations from python where we create 50.000 permutations
% % with alexander bloch 
% permutations = readmatrix([filepath '\DATA\spins_delta_permutations.csv']);
% 
% % extract the statistics of the linear mixed effect model for each feature
% subj=[];
% parcel=[];
% age=[];
% thick=[];
% y=[];
% for k1=1:350 %subjects
%     subj=[subj ones(1,200)*k1];
%     parcel=[parcel 1:200];
%     age = [age ones(1,200)*age_years(k1)];
%     thick=[thick thickness];
% end
% 
% 
% for rot  = 1:50000
%     perm=[];
%     for s = 1:350
%         perm =[perm permutations(:,rot)'];
%     end
% 
%     varnames={'subj','parcel','perm','age','thick'};
%     tbl=table(subj',parcel',perm',age',thick','VariableNames',varnames);
%     model = 'perm ~ thick*age + (1|subj)';
% 
%     try
%         lme=fitlme(tbl,model);
% 
%         %plotResiduals(lme,'fitted')
%         lm2(rot)=lme.Coefficients.tStat(3);
%     end
% end
% 
% 
% % find indices for 2.5% and 97.5% where the t-valzes should be significant
% idxmin5 = 0.025*50000;
% idxmin1 = 0.005*50000;
% idxmax5 = 0.975*50000;
% idxmax1 = 0.995*50000;
% 
% lm_sorted = sort(lm2,"ascend");
% lm_sorted(idxmax5) % 62.65
% lm_sorted(idxmax1) % 80.2
% lm_sorted(idxmin5) % -60.60
% lm_sorted(idxmin1) % -74.67
% 
% figure('Color','w', 'Position', [0 0 600 500])
% h = histogram(lm2,'NumBins',20)
% xlabel('LMEM t-values for permutated delta power map')
% ylabel('count')
% h.FaceColor = [0 0.5 0.5];
% 
% hold on 
% xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
% 
% xline(60.4,'r--','LineWidth',2,'Alpha',0.7) % Delta
% 
% box off
% 
% exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_delta_tvalues_thick.png'],'Resolution',500)
% 
% 
% %% THETA
% % use the permutations from python where we create 50.000 permutations
% % with alexander bloch 
% permutations = readmatrix([filepath '\DATA\spins_theta_permutations.csv']);
% 
% subj=[];
% parcel=[];
% age=[];
% thick=[];
% y=[];
% for k1=1:350 %subjects
%     subj=[subj ones(1,200)*k1];
%     parcel=[parcel 1:200];
%     age = [age ones(1,200)*age_years(k1)];
%     thick=[thick thickness];
% end
% 
% 
% for rot  = 1:50000
%     perm=[];
%     for s = 1:350
%         perm =[perm permutations(:,rot)'];
%     end
% 
%     varnames={'subj','parcel','perm','age','thick'};
%     tbl=table(subj',parcel',perm',age',thick','VariableNames',varnames);
%     model = 'perm ~ thick*age + (1|subj)';
% 
%     try
%         lme=fitlme(tbl,model);
% 
%         %plotResiduals(lme,'fitted')
%         lm3(rot)=lme.Coefficients.tStat(3);
%     end
% end
% 
% 
% % find indices for 2.5% and 97.5% where the t-valzes should be significant
% idxmin5 = 0.025*50000;
% idxmin1 = 0.005*50000;
% idxmax5 = 0.975*50000;
% idxmax1 = 0.995*50000;
% 
% lm_sorted = sort(lm3,"ascend");
% lm_sorted(idxmax5)
% lm_sorted(idxmax1) 
% lm_sorted(idxmin5)
% lm_sorted(idxmin1)
% 
% figure('Color','w', 'Position', [0 0 600 500])
% h = histogram(lm3,'NumBins',20)
% xlabel('LMEM t-values for permutated theta power map')
% ylabel('count')
% h.FaceColor = [0 0.5 0.5];
% 
% hold on 
% xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
% 
% xline(44.3,'r--','LineWidth',2,'Alpha',0.7) % t-value for theta power
% 
% box off
% 
% exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_theta_tvalues_thick.png'],'Resolution',500)
% 
% 
% 
% %% ALPHA
% % use the permutations from python where we create 50.000 permutations
% % with alexander bloch 
% permutations = readmatrix([filepath '\DATA\spins_alpha_permutations.csv']);
% 
% % extract the statistics of the linear mixed effect model for each feature
% subj=[];
% parcel=[];
% age=[];
% thick=[];
% y=[];
% for k1=1:350 %subjects
%     subj=[subj ones(1,200)*k1];
%     parcel=[parcel 1:200];
%     age = [age ones(1,200)*age_years(k1)];
%     thick=[thick thickness];
% end
% 
% 
% for rot  = 1:50000
%     perm=[];
%     for s = 1:350
%         perm =[perm permutations(:,rot)'];
%     end
% 
%     varnames={'subj','parcel','perm','age','thick'};
%     tbl=table(subj',parcel',perm',age',thick','VariableNames',varnames);
%     model = 'perm ~ thick*age + (1|subj)';
% 
%     try
%         lme=fitlme(tbl,model);
% 
%         %plotResiduals(lme,'fitted')
%         lm4(rot)=lme.Coefficients.tStat(3);
%     end
% end
% 
% 
% % find indices for 2.5% and 97.5% where the t-valzes should be significant
% idxmin5 = 0.025*50000;
% idxmin1 = 0.005*50000;
% idxmax5 = 0.975*50000;
% idxmax1 = 0.995*50000;
% 
% lm_sorted = sort(lm4,"ascend");
% lm_sorted(idxmax5)
% lm_sorted(idxmax1) 
% lm_sorted(idxmin5)
% lm_sorted(idxmin1)
% 
% figure('Color','w', 'Position', [0 0 600 500])
% h = histogram(lm4,'NumBins',20)
% xlabel('LMEM t-values for permutated alpha power map')
% ylabel('count')
% h.FaceColor = [0 0.5 0.5];
% 
% hold on 
% xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
% 
% xline(47.2,'r--','LineWidth',2,'Alpha',0.7) % t-value for alpha power
% 
% box off
% 
% exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_alpha_tvalues_thick.png'],'Resolution',500)
% 
% 
% 
% %% GAMMA
% % use the permutations from python where we create 50.000 permutations
% % with alexander bloch 
% permutations = readmatrix([filepath '\DATA\spins_gamma_permutations.csv']);
% 
% % extract the statistics of the linear mixed effect model for each feature
% subj=[];
% parcel=[];
% age=[];
% thick=[];
% y=[];
% for k1=1:350 %subjects
%     subj=[subj ones(1,200)*k1];
%     parcel=[parcel 1:200];
%     age = [age ones(1,200)*age_years(k1)];
%     thick=[thick thickness];
% end
% 
% 
% for rot  = 1:50000
%     perm=[];
%     for s = 1:350
%         perm =[perm permutations(:,rot)'];
%     end
% 
%     varnames={'subj','parcel','perm','age','thick'};
%     tbl=table(subj',parcel',perm',age',thick','VariableNames',varnames);
%     model = 'perm ~ thick*age + (1|subj)';
% 
%     try
%         lme=fitlme(tbl,model);
% 
%         %plotResiduals(lme,'fitted')
%         lm5(rot)=lme.Coefficients.tStat(3);
%     end
% end
% 
% % find indices for 2.5% and 97.5% where the t-valzes should be significant
% idxmin5 = 0.025*50000;
% idxmin1 = 0.005*50000;
% idxmax5 = 0.975*50000;
% idxmax1 = 0.995*50000;
% 
% lm_sorted = sort(lm5,"ascend");
% lm_sorted(idxmax5) 
% lm_sorted(idxmax1) 
% lm_sorted(idxmin5)
% lm_sorted(idxmin1)
% 
% figure('Color','w', 'Position', [0 0 600 500])
% h = histogram(lm5,'NumBins',20)
% xlabel('LMEM t-values for permutated gamma power map')
% ylabel('count')
% h.FaceColor = [0 0.5 0.5];
% 
% hold on 
% xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
% 
% xline(38.4,'r--','LineWidth',2,'Alpha',0.7) % t-value for gamma power
% 
% box off
% 
% exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_gamma_tvalues_thick.png'],'Resolution',500)
% 
% 
% 
% 
% 
% 





