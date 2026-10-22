%% create permutation for interaction effect

% we want to shuffle age in order to test our t-values form the interaction
% effect against a data set, where age should not have a systematic
% relationship with the change of time-series features across cortical
% hierarchy

filepath = pwd;

% load age
demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));

% load the parcellated map of cortical hierarchy (sydnor et al., 2021)
load([filepath '\DATA\hierarchy_parc.csv'])
 
% load mat_z file
load('G:\JF\HCTSA feature gradients\CamCan\big DATA files\mat_z.mat')

features = readtable([filepath '\DATA\features_z.csv']);


%% AC at lag 133 ms
idx_feature = find(strcmp("AC_40",table2cell(features(:,"Var2"))));
mat_feat = squeeze(mat_z(:,:,idx_feature));

% extract the statistics of the linear mixed effect model for each feature
subj=[];
parcel=[];
feat=[];
age=[];
hierarchy=[];

for k1=1:size(mat_z,1) %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    hierarchy=[hierarchy hierarchy_parc'];
end

for ishuffle=1:10000 % loop over features
    if (mod(ishuffle,200)==0)
        disp(ishuffle);
    end
    % shuffle age for each permutation round
    idx_shuffle = randsample(350,350);

    feat=[];
    age=[];
    for k1=1:size(mat_z,1) %subjects
        feat=[feat squeeze(mat_feat(k1,:))];
        age = [age ones(1,200)*age_years(idx_shuffle(k1))];
    end

    varnames={'subj','parcel','feat','age','hierarchy'};
    tbl=table(subj',parcel',feat',age',hierarchy','VariableNames',varnames);
    model = 'feat ~ hierarchy*age + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        allst_AC(ishuffle,:)=lme.Coefficients.tStat;
        pval_AC(ishuffle,:)=lme.Coefficients.pValue;
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(allst_AC(:,4),"ascend");
lm_sorted(idxmax5) 
lm_sorted(idxmax1) 
lm_sorted(idxmin5)
lm_sorted(idxmin1) 

% save allst_AC just in case
save([filepath '\DATA\surrogate_interarction_age_hierarchy_4AC40_results.mat'] ,"allst_AC")

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(allst_AC(:,4),'NumBins',20)
xlabel({'LMEM: t-values for AC (lag 133 ms)';'for the interaction effect of age and cortical hierarchy for shuffled age'})
ylabel('count')
h.FaceColor = [0 0.5 0.5];

hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)

xline(-42.7,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100

box off

exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_AC40_tValues_interactionAgeHierarchy_permAge.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
if sum(find(lm_sorted<-42.7)) == 0 
    p(1,1) = 1;
else
    p(1,1) =  (10000-max(find(lm_sorted<-42.7)))/10000;
end

%% plot AC40 across all ages
% 
% % plot
% parcellation = fetch_parcellation('fslr32k', 'schaefer', 200);
% [surf_lh, surf_rh] = fetch_template_surface('fslr32k', 'layer', 'inflated');
% load("vik.mat")
% 
% data2plot = mean(mat_feat,1);
% 
% vertex_data = nan(size(parcellation));
% for i = 1:200
%     vertex_data(parcellation == i) = data2plot(i);
% end
% obj = plot_hemispheres(vertex_data, {surf_lh, surf_rh});
% set(gcf,'units','points','position',[65,-65,1215,1215])
% colormap([0.7 0.7 0.7;vik])
% saveas(gcf,'x.png')
% I1 = imread('x.png');
% T1 = imcrop(I1, [180, 60, 555, 210]); % only left hemisphere
% Imerged = imtile({T1},  'BorderSize', 4, 'BackgroundColor', 'white','GridSize', [1,1]);
% imwrite(Imerged,[filepath '\PLOTS\SUPPL_topo_AC40.png']);



%% AMI40
load([filepath '\DATA\ami_1_40.mat'])
ami = squeeze(ami_out(:,:,40));
%zscore
ami_z= zscore(ami,[],2);

% extract the statistics of the linear mixed effect model for each feature
subj=[];
parcel=[];
feat=[];
age=[];
hierarchy=[];

for k1=1:size(mat_z,1) %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    hierarchy=[hierarchy hierarchy_parc'];
end

for ishuffle=1:10000 % loop over shuffels
    if (mod(ishuffle,200)==0)
        disp(ishuffle);
    end
    % shuffle age for each permutation round
    idx_shuffle = randsample(350,350);

    feat=[];
    age=[];
    for k1=1:size(mat_z,1) %subjects
        feat=[feat ami_z(k1,:)];
        age = [age ones(1,200)*age_years(idx_shuffle(k1))];
    end

    varnames={'subj','parcel','feat','age','hierarchy'};
    tbl=table(subj',parcel',feat',age',hierarchy','VariableNames',varnames);
    model = 'feat ~ hierarchy*age + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        allst_AMI(ishuffle,:)=lme.Coefficients.tStat;
        pval_AMI(ishuffle,:)=lme.Coefficients.pValue;
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(allst_AMI(:,4),"ascend");
lm_sorted(idxmax5) % 5.20
lm_sorted(idxmax1) % 6.87
lm_sorted(idxmin5) % -5.21
lm_sorted(idxmin1) % -6.71

% save allst_AC just in case
save([filepath '\DATA\surrogate_interarction_age_hierarchy_4AMI40_results.mat'],"allst_AMI")

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(allst_AMI(:,4),'NumBins',20)
xlabel({'LMEM: t-values for AMI (lag 133 ms)';'for the interaction effect of age and cortical hierarchy for shuffled age'})
ylabel('count')
h.FaceColor = [0 0.5 0.5];

hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)

xline(27.65,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100

box off

exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_AMI40_tValues_interactionAgeHierarchy_permAge.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
if sum(find(lm_sorted<27.65)) == 0 
    p(1,2) = 1;
else
    p(1,2) =  (10000-max(find(lm_sorted<27.65)))/10000;
end


% % plot
% parcellation = fetch_parcellation('fslr32k', 'schaefer', 200);
% [surf_lh, surf_rh] = fetch_template_surface('fslr32k', 'layer', 'inflated');
% load("vik.mat")
% 
% data2plot = mean(ami_z,1)';
% 
% vertex_data = nan(size(parcellation));
% for i = 1:200
%     vertex_data(parcellation == i) = data2plot(i);
% end
% obj = plot_hemispheres(vertex_data, {surf_lh, surf_rh});
% set(gcf,'units','points','position',[65,-65,1215,1215])
% colormap([0.7 0.7 0.7;vik])
% saveas(gcf,'x.png')
% I1 = imread('x.png');
% T1 = imcrop(I1, [180, 60, 555, 210]); % only left hemisphere
% Imerged = imtile({T1},  'BorderSize', 4, 'BackgroundColor', 'white','GridSize', [1,1]);
% imwrite(Imerged,[filepath '\PLOTS\SUPPL_topo_AMI40.png']);
% 
% 

%% DELTA 
% loading delta
load([filepath '\DATA\03_specData_freqbands_norm'],"delta_norm");
delta_z = zscore(delta_norm,[],2);

% extract the statistics of the linear mixed effect model for each feature
subj=[];
parcel=[];
feat=[];
age=[];
hierarchy=[];

for k1=1:350 %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    hierarchy=[hierarchy hierarchy_parc'];
end

for ishuffle=1:10000 % loop over features
    if (mod(ishuffle,200)==0)
        disp(ishuffle);
    end
    % shuffle age for each permutation round
    idx_shuffle = randsample(350,350);

    feat=[];
    age=[];
    for k1=1:350 %subjects
        feat=[feat delta_z(k1,:)];
        age = [age ones(1,200)*age_years(idx_shuffle(k1))];
    end

    varnames={'subj','parcel','feat','age','hierarchy'};
    tbl=table(subj',parcel',feat',age',hierarchy','VariableNames',varnames);
    model = 'feat ~ hierarchy*age + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        allst_delta(ishuffle,:)=lme.Coefficients.tStat;
        pval_delta(ishuffle,:)=lme.Coefficients.pValue;
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(allst_delta(:,4),"ascend");
lm_sorted(idxmax5)
lm_sorted(idxmax1) 
lm_sorted(idxmin5) 
lm_sorted(idxmin1) 

% save allst_AC just in case
save([filepath '\DATA\surrogate_interarction_age_hierarchy_4delta_results.mat'],"allst_delta")

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(allst_delta(:,4),'NumBins',20)
xlabel({'LMEM: t-values for delta power';'for the interaction effect of age and cortical hierarchy for shuffled age'})
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(-17.62,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_delta_tValues_interactionAgeHierarchy_permAge.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
if sum(find(lm_sorted<-17.62)) == 0 
    p(1,3) = 1;
else
    p(1,3) =  (10000-max(find(lm_sorted<-17.62)))/10000;
end



%% THETA 
% loading theta
load([filepath '\DATA\03_specData_freqbands_norm'],"theta_norm");
theta_z = zscore(theta_norm,[],2);

for ishuffle=1:10000 % loop over features
    if (mod(ishuffle,200)==0)
        disp(ishuffle);
    end
    % shuffle age for each permutation round
    idx_shuffle = randsample(350,350);

    feat=[];
    age=[];
    for k1=1:350 %subjects
        feat=[feat theta_z(k1,:)];
        age = [age ones(1,200)*age_years(idx_shuffle(k1))];
    end

    varnames={'subj','parcel','feat','age','hierarchy'};
    tbl=table(subj',parcel',feat',age',hierarchy','VariableNames',varnames);
    model = 'feat ~ hierarchy*age + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        allst_theta(ishuffle,:)=lme.Coefficients.tStat;
        pval_theta(ishuffle,:)=lme.Coefficients.pValue;
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(allst_theta(:,4),"ascend");
lm_sorted(idxmax5) 
lm_sorted(idxmax1) 
lm_sorted(idxmin5)
lm_sorted(idxmin1) 

% save allst_AC just in case
save([filepath '\DATA\surrogate_interarction_age_hierarchy_theta_results.mat'],"allst_theta")

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(allst_theta(:,4),'NumBins',20)
xlabel({'LMEM: t-values for theta power';'for the interaction effect of age and cortical hierarchy for shuffled age'})
ylabel('count')
h.FaceColor = [0 0.5 0.5];

hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)

xline(-18.63,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100

box off

exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_theta_tValues_interactionAgeHierarchy_permAge.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
if sum(find(lm_sorted<-18.63)) == 0 
    p(1,4) = 1;
else
    p(1,4) =  (10000-max(find(lm_sorted<-18.63)))/10000;
end

%% ALPHA 
% loading theta
load([filepath '\DATA\03_specData_freqbands_norm'],"alpha_norm");
alpha_z = zscore(alpha_norm,[],2);

% extract the statistics of the linear mixed effect model for each feature
subj=[];
parcel=[];
feat=[];
age=[];
hierarchy=[];

for k1=1:size(mat_z,1) %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    hierarchy=[hierarchy hierarchy_parc'];
end

for ishuffle=1:10000 % loop over features
    if (mod(ishuffle,200)==0)
        disp(ishuffle);
    end
    % shuffle age for each permutation round
    idx_shuffle = randsample(350,350);

    feat=[];
    age=[];
    for k1=1:350 %subjects
        feat=[feat alpha_z(k1,:)];
        age = [age ones(1,200)*age_years(idx_shuffle(k1))];
    end

    varnames={'subj','parcel','feat','age','hierarchy'};
    tbl=table(subj',parcel',feat',age',hierarchy','VariableNames',varnames);
    model = 'feat ~ hierarchy*age + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        allst_alpha(ishuffle,:)=lme.Coefficients.tStat;
        pval_alpha(ishuffle,:)=lme.Coefficients.pValue;
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(allst_alpha(:,4),"ascend");
lm_sorted(idxmax5) 
lm_sorted(idxmax1)
lm_sorted(idxmin5) 
lm_sorted(idxmin1) % 

% save allst_AC just in case
save([filepath '\DATA\surrogate_interarction_age_hierarchy_4alpha_results.mat'],"allst_alpha")

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(allst_alpha(:,4),'NumBins',20)
xlabel({'LMEM: t-values for alpha power';'for the interaction effect of age and cortical hierarchy for shuffled age'})
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(38.05,'r--','LineWidth',2,'Alpha',0.7) 
box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_alpha_tValues_interactionAgeHierarchy_permAge.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
if sum(find(lm_sorted<38.05)) == 0 
    p(1,5) = 1;
else
    p(1,5) =  (10000-max(find(lm_sorted<38.05)))/10000;
end




%% BETA 
% loading theta
load([filepath '\DATA\03_specData_freqbands_norm'],"beta_norm");
beta_z = zscore(beta_norm,[],2);

% extract the statistics of the linear mixed effect model for each feature

feat=[];
age=[];

for ishuffle=1:10000 % loop over features
    if (mod(ishuffle,200)==0)
        disp(ishuffle);
    end
    % shuffle age for each permutation round
    idx_shuffle = randsample(350,350);

    feat=[];
    age=[];
    for k1=1:350 %subjects
        feat=[feat beta_z(k1,:)];
        age = [age ones(1,200)*age_years(idx_shuffle(k1))];
    end

    varnames={'subj','parcel','feat','age','hierarchy'};
    tbl=table(subj',parcel',feat',age',hierarchy','VariableNames',varnames);
    model = 'feat ~ hierarchy*age + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        allst_beta(ishuffle,:)=lme.Coefficients.tStat;
        pval_beta(ishuffle,:)=lme.Coefficients.pValue;
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(allst_beta(:,4),"ascend");
lm_sorted(idxmax5) 
lm_sorted(idxmax1) 
lm_sorted(idxmin5) 
lm_sorted(idxmin1)

% save allst_AC just in case
save([filepath '\DATA\surrogate_interarction_age_hierarchy_beta_results.mat'],"allst_beta")

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(allst_beta(:,4),'NumBins',20)
xlabel({'LMEM: t-values for beta power';'for the interaction effect of age and cortical hierarchy for shuffled age'})
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(-7.72,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_beta_tValues_interactionAgeHierarchy_permAge.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
p(1,6) =  (10000-max(find(lm_sorted<-7.72)))/10000;




%% GAMMA 
% loading theta
load([filepath '\DATA\03_specData_freqbands_norm'],"gamma_norm");
gamma_z = zscore(gamma_norm,[],2);

% extract the statistics of the linear mixed effect model for each feature
subj=[];
parcel=[];
feat=[];
age=[];
hierarchy=[];

for k1=1:350 %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    hierarchy=[hierarchy hierarchy_parc'];
end

for ishuffle=1:10000 % loop over features
    if (mod(ishuffle,200)==0)
        disp(ishuffle);
    end
    % shuffle age for each permutation round
    idx_shuffle = randsample(350,350);

    feat=[];
    age=[];
    for k1=1:350 %subjects
        feat=[feat gamma_z(k1,:)];
        age = [age ones(1,200)*age_years(idx_shuffle(k1))];
    end

    varnames={'subj','parcel','feat','age','hierarchy'};
    tbl=table(subj',parcel',feat',age',hierarchy','VariableNames',varnames);
    model = 'feat ~ hierarchy*age + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        allst_gamma(ishuffle,:)=lme.Coefficients.tStat;
        pval_gamma(ishuffle,:)=lme.Coefficients.pValue;
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(allst_gamma(:,4),"ascend");
lm_sorted(idxmax5)
lm_sorted(idxmax1)
lm_sorted(idxmin5)
lm_sorted(idxmin1)

% save allst_AC just in case
save([filepath '\DATA\surrogate_interarction_age_hierarchy_4gamma_results.mat'],"allst_gamma")

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(allst_gamma(:,4),'NumBins',20)
xlabel({'LMEM: t-values for gamma power';'for the interaction effect of age and cortical hierarchy for shuffled age'})
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(-18.16,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_gamma_tValues_interactionAgeHierarchy_permAge.png'],'Resolution',500)


% which p-value does the t-vual correspond to?
if sum(find(lm_sorted<-18.16)) == 0 
    p(1,7) = 1;
else
    p(1,7) =  (10000-max(find(lm_sorted<-18.16)))/10000;
end



%% CENTRE OF GRAVITY
% loading theta
load([filepath '\DATA\peak_freq_cog.mat'])
cog = zscore(peak_freq_cog,[],2);

% extract the statistics of the linear mixed effect model for each feature
subj=[];
parcel=[];
feat=[];
age=[];
hierarchy=[];

for k1=1:350 %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    hierarchy=[hierarchy hierarchy_parc'];
end

for ishuffle=1:10000 % loop over features
    if (mod(ishuffle,200)==0)
        disp(ishuffle);
    end
    % shuffle age for each permutation round
    idx_shuffle = randsample(350,350);

    feat=[];
    age=[];
    for k1=1:350 %subjects
        feat=[feat cog(k1,:)];
        age = [age ones(1,200)*age_years(idx_shuffle(k1))];
    end

    varnames={'subj','parcel','feat','age','hierarchy'};
    tbl=table(subj',parcel',feat',age',hierarchy','VariableNames',varnames);
    model = 'feat ~ hierarchy*age + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        allst_cog(ishuffle,:)=lme.Coefficients.tStat;
        pval_cog(ishuffle,:)=lme.Coefficients.pValue;
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(allst_cog(:,4),"ascend");
lm_sorted(idxmax5) 
lm_sorted(idxmax1) 
lm_sorted(idxmin5) 
lm_sorted(idxmin1) 

% save allst_AC just in case
save([filepath '\DATA\surrogate_interarction_age_hierarchy_4cog_results.mat'],"allst_cog")

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(allst_cog(:,4),'NumBins',20)
xlabel({'LMEM: t-values for peak frequency';'for the interaction effect of age and cortical hierarchy for shuffled age'})
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(18.96,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_cog_tValues_interactionAgeHierarchy_permAge.png'],'Resolution',500)


% which p-value does the t-vual correspond to?
if sum(find(lm_sorted<18.96)) == 0 
    p(1,8) = 1;
else
    p(1,8) =  (10000-max(find(lm_sorted<18.96)))/10000;
end



%% FOOOF
% loading theta
load([filepath '\DATA\correctedFooofParam_fooofWknee_z.mat']);
fooof = exp_z;

% extract the statistics of the linear mixed effect model for each feature
subj=[];
parcel=[];
feat=[];
age=[];
hierarchy=[];

for k1=1:size(mat_z,1) %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    hierarchy=[hierarchy hierarchy_parc'];
end

for ishuffle=1:10000 % loop over features
    if (mod(ishuffle,200)==0)
        disp(ishuffle);
    end
    % shuffle age for each permutation round
    idx_shuffle = randsample(350,350);

    feat=[];
    age=[];
    for k1=1:size(mat_z,1) %subjects
        feat=[feat fooof(k1,:)];
        age = [age ones(1,200)*age_years(idx_shuffle(k1))];
    end

    varnames={'subj','parcel','feat','age','hierarchy'};
    tbl=table(subj',parcel',feat',age',hierarchy','VariableNames',varnames);
    model = 'feat ~ hierarchy*age + (1|subj)';
    try
        lme=fitlme(tbl,model);
        allst_fooof(ishuffle,:)=lme.Coefficients.tStat;
        pval_fooof(ishuffle,:)=lme.Coefficients.pValue;
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(allst_fooof(:,4),"ascend");
lm_sorted(idxmax5) % 4.46
lm_sorted(idxmax1) % 5.77
lm_sorted(idxmin5) % -4.55
lm_sorted(idxmin1) % -5.91

% save allst_AC just in case
save([filepath '\DATA\surrogate_interarction_age_hierarchy_4fooof_results.mat'],"allst_fooof")

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(allst_fooof(:,4),'NumBins',20)
xlabel({'LMEM: t-values for 1/f component';'for the interaction effect of age and cortical hierarchy for shuffled age'})
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(25.98,'r--','LineWidth',2,'Alpha',0.7)
box off

exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_fooof_tValues_interactionAgeHierarchy_permAge.png'],'Resolution',500)


% which p-value does the t-vual correspond to?
if sum(find(lm_sorted<25.98)) == 0 
    p(1,9) = 1;
else
    p(1,9) =  (10000-max(find(lm_sorted<25.98)))/10000;
end








%%%%%%%%% NOW DO THE SAME WITH CORTICAL HIERARCHY %%%%%%%%%%%%%%%%%%%%

% loading cortical thickness
load([filepath '\DATA\thickness_smoothed_12_schaefer200_reordered.mat']);
thickness(306,:) = [];

% delete participant 306 for all variables --> huge cortical thickness
% outlier
age_years(306) = [];

% extract the statistics of the linear mixed effect model for each feature
subj=[];
parcel=[];
thick=[];
y=[];
for k1=1:349 %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    thick=[thick thickness(k1,:)];
end


%% AUTOCORRELATION 40
idx_feature = find(strcmp("AC_40",table2cell(features(:,"Var2"))));
mat_feat = squeeze(mat_z(:,:,idx_feature));
mat_feat(306,:) = [];

% extract the statistics of the linear mixed effect model for each feature
feat=[];
age=[];

for ishuffle=1:10000 % loop over features
    if (mod(ishuffle,200)==0)
        disp(ishuffle);
    end
    % shuffle age for each permutation round
    idx_shuffle = randsample(349,349);

    feat=[];
    age=[];
    for k1=1:349 %subjects
        feat=[feat squeeze(mat_feat(k1,:))];
        age = [age ones(1,200)*age_years(idx_shuffle(k1))];
    end

    varnames={'subj','parcel','feat','age','thick'};
    tbl=table(subj',parcel',feat',age',thick','VariableNames',varnames);
    model = 'feat ~ thick*age + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        allst_AC(ishuffle,:)=lme.Coefficients.tStat;
        pval_AC(ishuffle,:)=lme.Coefficients.pValue;
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(allst_AC(:,4),"ascend");
lm_sorted(idxmax5) 
lm_sorted(idxmax1) 
lm_sorted(idxmin5) 
lm_sorted(idxmin1) 

% save allst_AC just in case
save([filepath '\DATA\surrogate_interarction_age_thickness_4AC40_results.mat'] ,"allst_AC")

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(allst_AC(:,4),'NumBins',20)
xlabel({'LMEM: t-values for AC (lag 133 ms)';'for the interaction effect of age and cortical thickness for shuffled age'})
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(-19.25,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_AC40_tValues_interactionAgeThick_permAge.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
if sum(find(lm_sorted<-19.25)) == 0 
    p(2,1) = 1;
else
    p(2,1) =  (10000-max(find(lm_sorted<-19.25)))/10000;
end


%% AUTO MUTUAL INFORMATION 40
load([filepath '\DATA\ami_1_40.mat'])
ami = squeeze(ami_out(:,:,40));
%zscore
ami_z= zscore(ami,[],2); ami_z(306,:);

% extract the statistics of the linear mixed effect model for each feature
feat=[];
age=[];

for ishuffle=1:10000 % loop over features
    if (mod(ishuffle,200)==0)
        disp(ishuffle);
    end
    % shuffle age for each permutation round
    idx_shuffle = randsample(349,349);

    feat=[];
    age=[];
    for k1=1:349 %subjects
        feat=[feat squeeze(ami_z(k1,:))];
        age = [age ones(1,200)*age_years(idx_shuffle(k1))];
    end

    varnames={'subj','parcel','feat','age','thick'};
    tbl=table(subj',parcel',feat',age',thick','VariableNames',varnames);
    model = 'feat ~ thick*age + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        allst_AMI(ishuffle,:)=lme.Coefficients.tStat;
        pval_AMI(ishuffle,:)=lme.Coefficients.pValue;
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(allst_AMI(:,4),"ascend");
lm_sorted(idxmax5) % 9.15
lm_sorted(idxmax1) % 12.05
lm_sorted(idxmin5) % -9.06
lm_sorted(idxmin1) % -11.68

% save allst_AC just in case
save([filepath '\DATA\surrogate_interarction_age_thickness_4MI40_results.mat'] ,"allst_AMI")

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(allst_AMI(:,4),'NumBins',20)
xlabel({'LMEM: t-values for AMI (lag 133 ms)';'for the interaction effect of age and cortical thickness for shuffled age'})
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(26.33,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_AMI40_tValues_interactionAgeThick_permAge.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
if sum(find(lm_sorted<26,33)) == 0 
    p(2,2) = 1;
else
    p(2,2) =  (10000-max(find(lm_sorted<26.33)))/10000;
end


%% DELTA 

% loading delta
load([filepath '\DATA\03_specData_freqbands_norm'],"delta_norm");
delta_z = zscore(delta_norm,[],2); delta_z(306,:)=[];

% extract the statistics of the linear mixed effect model for each feature
feat=[];
age=[];

for ishuffle=1:10000 % loop over features
    if (mod(ishuffle,200)==0)
        disp(ishuffle);
    end
    % shuffle age for each permutation round
    idx_shuffle = randsample(349,349);

    feat=[];
    age=[];
    for k1=1:349 %subjects
        feat=[feat squeeze(delta_z(k1,:))];
        age = [age ones(1,200)*age_years(idx_shuffle(k1))];
    end

    varnames={'subj','parcel','feat','age','thick'};
    tbl=table(subj',parcel',feat',age',thick','VariableNames',varnames);
    model = 'feat ~ thick*age + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        allst(ishuffle,:)=lme.Coefficients.tStat;
        pval(ishuffle,:)=lme.Coefficients.pValue;
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(allst(:,4),"ascend");
lm_sorted(idxmax5) % 9.15
lm_sorted(idxmax1) % 12.05
lm_sorted(idxmin5) % -9.06
lm_sorted(idxmin1) % -11.68

% save allst just in case
save([filepath '\DATA\surrogate_interarction_age_thickness_delta_results.mat'] ,"allst")



figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm_sorted,'NumBins',20)
xlabel({'LMEM: t-values for delta power';'for the interaction effect of age and cortical thickness for shuffled age'})
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(-22.92,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off

exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_delta_tValues_interactionAgeThick_permAge.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
if sum(find(lm_sorted<-22.92)) == 0 
    p(2,3) = 1;
else
    p(2,3) =  (10000-max(find(lm_sorted<-22.92)))/10000;
end


%% THETA 

load([filepath '\DATA\03_specData_freqbands_norm'],"theta_norm");
theta_z = zscore(theta_norm,[],2); theta_z(306,:)=[];

% extract the statistics of the linear mixed effect model for each feature
feat=[];
age=[];

for ishuffle=1:10000 % loop over features
    if (mod(ishuffle,200)==0)
        disp(ishuffle);
    end
    % shuffle age for each permutation round
    idx_shuffle = randsample(349,349);

    feat=[];
    age=[];
    for k1=1:349 %subjects
        feat=[feat squeeze(theta_z(k1,:))];
        age = [age ones(1,200)*age_years(idx_shuffle(k1))];
    end

    varnames={'subj','parcel','feat','age','thick'};
    tbl=table(subj',parcel',feat',age',thick','VariableNames',varnames);
    model = 'feat ~ thick*age + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        allst(ishuffle,:)=lme.Coefficients.tStat;
        pval(ishuffle,:)=lme.Coefficients.pValue;
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(allst(:,4),"ascend");
lm_sorted(idxmax5) 
lm_sorted(idxmax1) 
lm_sorted(idxmin5) 
lm_sorted(idxmin1) 

% save allst just in case
save([filepath '\DATA\surrogate_interarction_age_thickness_theta_results.mat'] ,"allst")


figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm_sorted,'NumBins',20)
xlabel({'LMEM: t-values for theta power';'for the interaction effect of age and cortical thickness for shuffled age'})
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(-14.23,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_theta_tValues_interactionAgeThick_permAge.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
if sum(find(lm_sorted<-14.23)) == 0 
    p(2,4) = 1;
else
    p(2,4) =  (10000-max(find(lm_sorted<-14.23)))/10000;
end




%% ALPHA 

load([filepath '\DATA\03_specData_freqbands_norm'],"alpha_norm");
alpha_z = zscore(alpha_norm,[],2);alpha_z(306,:)=[];


for ishuffle=1:10000 % loop over features
    if (mod(ishuffle,200)==0)
        disp(ishuffle);
    end
    % shuffle age for each permutation round
    idx_shuffle = randsample(349,349);

    feat=[];
    age=[];
    for k1=1:349 %subjects
        feat=[feat alpha_z(k1,:)];
        age = [age ones(1,200)*age_years(idx_shuffle(k1))];
    end

    varnames={'subj','parcel','feat','age','thick'};
    tbl=table(subj',parcel',feat',age',thick','VariableNames',varnames);
    model = 'feat ~ thick*age + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        allst(ishuffle,:)=lme.Coefficients.tStat;
        pval(ishuffle,:)=lme.Coefficients.pValue;
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(allst(:,4),"ascend");
lm_sorted(idxmax5) 
lm_sorted(idxmax1) 
lm_sorted(idxmin5) 
lm_sorted(idxmin1) 

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm_sorted,'NumBins',20)
xlabel({'LMEM: t-values for alpha power';'for the interaction effect of age and cortical thickness for shuffled age'})
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(33.3,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_alpha_tValues_interactionAgeThick_permAge.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
if sum(find(lm_sorted<33.3)) == 0 
    p(2,5) = 1;
else
    p(2,5) =  (10000-max(find(lm_sorted<33.3)))/10000;
end




%% BETA

load([filepath '\DATA\03_specData_freqbands_norm'],"beta_norm");
beta_z = zscore(alpha_norm,[],2);beta_z(306,:)=[];

for ishuffle=1:10000 % loop over features
    if (mod(ishuffle,200)==0)
        disp(ishuffle);
    end
    % shuffle age for each permutation round
    idx_shuffle = randsample(349,349);

    feat=[];
    age=[];
    for k1=1:349 %subjects
        feat=[feat beta_z(k1,:)];
        age = [age ones(1,200)*age_years(idx_shuffle(k1))];
    end

    varnames={'subj','parcel','feat','age','thick'};
    tbl=table(subj',parcel',feat',age',thick','VariableNames',varnames);
    model = 'feat ~ thick*age + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        allst(ishuffle,:)=lme.Coefficients.tStat;
        pval(ishuffle,:)=lme.Coefficients.pValue;
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(allst(:,4),"ascend");
lm_sorted(idxmax5) 
lm_sorted(idxmax1) 
lm_sorted(idxmin5) 
lm_sorted(idxmin1) 

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm_sorted,'NumBins',20)
xlabel({'LMEM: t-values for beta power';'for the interaction effect of age and cortical thickness for shuffled age'})
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(-0.73,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off

exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_beta_tValues_interactionAgethick_permAge.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
if sum(find(lm_sorted<-0.73)) == 0 
    p(2,6) = 1;
else
    p(2,6) =  (10000-max(find(lm_sorted<-0.73)))/10000;
end




%% GAMMA 

load([filepath '\DATA\03_specData_freqbands_norm'],"gamma_norm");
gamma_z = zscore(gamma_norm,[],2);gamma_z(306,:)=[];

for ishuffle=1:10000 % loop over features
    if (mod(ishuffle,200)==0)
        disp(ishuffle);
    end
    % shuffle age for each permutation round
    idx_shuffle = randsample(349,349);

    feat=[];
    age=[];
    for k1=1:349 %subjects
        feat=[feat gamma_z(k1,:)];
        age = [age ones(1,200)*age_years(idx_shuffle(k1))];
    end

    varnames={'subj','parcel','feat','age','thick'};
    tbl=table(subj',parcel',feat',age',thick','VariableNames',varnames);
    model = 'feat ~ thick*age + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        allst(ishuffle,:)=lme.Coefficients.tStat;
        pval(ishuffle,:)=lme.Coefficients.pValue;
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(allst(:,4),"ascend");
lm_sorted(idxmax5) 
lm_sorted(idxmax1) 
lm_sorted(idxmin5) 
lm_sorted(idxmin1) 

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm_sorted,'NumBins',20)
xlabel({'LMEM: t-values for gamma power';'for the interaction effect of age and cortical thickness for shuffled age'})
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(-17.02,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_gamma_tValues_interactionAgethick_permAge.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
if sum(find(lm_sorted<-17.02)) == 0 
    p(2,7) = 1;
else
    p(2,7) =  (10000-max(find(lm_sorted<-17.02)))/10000;
end



%%  cog 

% loading theta
load([filepath '\DATA\foof_variables.mat'])
cog = zscore(peak_freq_cog,[],2);cog(306,:)=[];

for ishuffle=1:10000 % loop over features
    if (mod(ishuffle,200)==0)
        disp(ishuffle);
    end
    % shuffle age for each permutation round
    idx_shuffle = randsample(349,349);

    feat=[];
    age=[];
    for k1=1:349 %subjects
        feat=[feat cog(k1,:)];
        age = [age ones(1,200)*age_years(idx_shuffle(k1))];
    end

    varnames={'subj','parcel','feat','age','thick'};
    tbl=table(subj',parcel',feat',age',thick','VariableNames',varnames);
    model = 'feat ~ thick*age + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        allst(ishuffle,:)=lme.Coefficients.tStat;
        pval(ishuffle,:)=lme.Coefficients.pValue;
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(allst(:,4),"ascend");
lm_sorted(idxmax5) 
lm_sorted(idxmax1) 
lm_sorted(idxmin5) 
lm_sorted(idxmin1) 

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm_sorted,'NumBins',20)
xlabel({'LMEM: t-values for peak frequency';'for the interaction effect of age and cortical thickness for shuffled age'})
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(5.17,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_cog_tValues_interactionAgethick_permAge.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
if sum(find(lm_sorted<5.17)) == 0 
    p(2,8) = 1;
else
    p(2,8) =  (10000-max(find(lm_sorted<5.17)))/10000;
end



%% FOOOF
% loading theta
load([filepath '\DATA\correctedFooofParam_fooofWknee_z.mat']);
fooof = exp_z; fooof(306,:) = [];

for ishuffle=1:10000 % loop over features
    if (mod(ishuffle,200)==0)
        disp(ishuffle);
    end
    % shuffle age for each permutation round
    idx_shuffle = randsample(349,349);

    feat=[];
    age=[];
    for k1=1:349 %subjects
        feat=[feat fooof(k1,:)];
        age = [age ones(1,200)*age_years(idx_shuffle(k1))];
    end

    varnames={'subj','parcel','feat','age','thick'};
    tbl=table(subj',parcel',feat',age',thick','VariableNames',varnames);
    model = 'feat ~ thick*age + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        allst(ishuffle,:)=lme.Coefficients.tStat;
        pval(ishuffle,:)=lme.Coefficients.pValue;
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*10000;
idxmin1 = 0.005*10000;
idxmax5 = 0.975*10000;
idxmax1 = 0.995*10000;

lm_sorted = sort(allst(:,4),"ascend");
lm_sorted(idxmax5) % 4.46
lm_sorted(idxmax1) % 5.77
lm_sorted(idxmin5) % -4.55
lm_sorted(idxmin1) % -5.91

% save allst_AC just in case
save([filepath '\DATA\surrogate_interarction_age_hierarchy_4fooof_results.mat'],"allst_fooof")

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(allst_fooof(:,4),'NumBins',20)
xlabel({'LMEM: t-values for 1/f component';'for the interaction effect of age and cortical thickness for shuffled age'})
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(20.81,'r--','LineWidth',2,'Alpha',0.7)
box off

exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_fooof_tValues_interactionAgeThick_permAge.png'],'Resolution',500)


% which p-value does the t-vual correspond to?
if sum(find(lm_sorted<20.81)) == 0 
    p(2,9) = 1;
else
    p(2,9) =  (10000-max(find(lm_sorted<20.81)))/10000;
end
