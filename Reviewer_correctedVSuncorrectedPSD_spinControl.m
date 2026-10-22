% CONTROL ANALYSIS 
% control analysis for LMEM: spinning corticla hierarchy or corticla
% thickness to get a distribution of expectable t-values to test original t-value against
clear all
close all

filepath = pwd;

% load the parcellated map of cortical hierarchy (sydnor et al., 2021)
load([filepath '\DATA\hierarchy_parc.csv']);

% load demogrpahics in order to catch the age in years
demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));      

% load Cortcial hierarchy permutations (50.000 according to alexander bloch)
permutations = readmatrix([filepath '\DATA\spins_hierarchy_permutations.csv']);

load([filepath '\DATA\correctedPowSpec_fooofWknee_norm_z.mat'],'corrPowSpec');
load([filepath '\DATA\correctedFooofParam_fooofWknee_z.mat'],'exp_z',"offs_z");

%% DELTA
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
    feat = [feat squeeze(corrPowSpec(1,k1,:))'];
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
xlabel('LMEM t-values for permutated delta map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)

xline(76.96,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100

box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_corrdelta_tvalues.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
p(1,1) =  (10000-max(find(lm_sorted<76.96)))/10000


%% THETA
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
    feat = [feat squeeze(corrPowSpec(2,k1,:))'];
end

for rot  = 1:10000
    if ~mod(rot,1000)
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
lm_sorted(idxmax5) % 
lm_sorted(idxmax1) % 
lm_sorted(idxmin5) % 
lm_sorted(idxmin1) % 

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm,'NumBins',20)
xlabel('LMEM t-values for permutated delta map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(58.08,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_corrtheta_tvalues.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
p(1,2) =  (10000-max(find(lm_sorted<58.08)))/10000


%% ALPHA
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
    feat = [feat squeeze(corrPowSpec(3,k1,:))']; 
end

for rot  = 1:10000
    if ~mod(rot,1000)
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
lm_sorted(idxmax5) % 
lm_sorted(idxmax1) % 
lm_sorted(idxmin5) % 
lm_sorted(idxmin1) % 

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm,'NumBins',20)
xlabel('LMEM t-values for permutated delta map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)

xline(-63.68,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100

box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_corralpha_tvalues.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
p(1,3) =  (10000-max(find(lm_sorted<-63.68)))/10000


%% BETA
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
    feat = [feat squeeze(corrPowSpec(4,k1,:))']; 
end

for rot  = 1:10000
    if ~mod(rot,1000)
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
lm_sorted(idxmax5) % 
lm_sorted(idxmax1) % 
lm_sorted(idxmin5) % 
lm_sorted(idxmin1) % 

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm,'NumBins',20)
xlabel('LMEM t-values for permutated delta map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(-5.08,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_corrbeta_tvalues.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
p(1,4) =  (10000-max(find(lm_sorted<-5.08)))/10000


%% GAMMA
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
    feat = [feat squeeze(corrPowSpec(5,k1,:))']; 
end

for rot  = 1:10000
    if ~mod(rot,1000)
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
lm_sorted(idxmax5) % 
lm_sorted(idxmax1) % 
lm_sorted(idxmin5) % 
lm_sorted(idxmin1) % 

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm,'NumBins',20)
xlabel('LMEM t-values for permutated delta map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(64.4,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_corrgamma_tvalues.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
p(1,5) =  (10000-max(find(lm_sorted<64.4)))/10000


%% FOOOF EXPONENT
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
    feat = [feat exp_z(k1,:)]; 
end

for rot  = 1:10000
    if ~mod(rot,1000)
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
lm_sorted(idxmax5) % 
lm_sorted(idxmax1) % 
lm_sorted(idxmin5) % 
lm_sorted(idxmin1) % 

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm,'NumBins',20)
xlabel('LMEM t-values for permutated aperiodic exponent map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(-14.18,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_fooofexp_tvalues.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
p(1,6) =  (10000-max(find(lm_sorted<-14.18)))/10000


%% FOOOF OFFSET
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
    feat = [feat offs_z(k1,:)]; 
end

for rot  = 1:10000
    if ~mod(rot,1000)
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
lm_sorted(idxmax5) % 
lm_sorted(idxmax1) % 
lm_sorted(idxmin5) % 
lm_sorted(idxmin1) % 

figure('Color','w', 'Position', [0 0 600 500])
h = histogram(lm,'NumBins',20)
xlabel('LMEM t-values for permutated delta map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(-66.91,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_fooofoffs_tvalues.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
p(1,7) =  (10000-max(find(lm_sorted<-66.91)))/10000




%% ====================================================================
% CORTICAL THHICKNESS
% do thesame with cortical thickness
spin_index = readmatrix([filepath, '\DATA\spins_index_permutations.csv']);

% loading cortical thickness
load([filepath '\DATA\thickness_smoothed_12_schaefer200_reordered.mat']);
thickness(306,:) = [];

demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));     
age_years(306,:) = [];   

corrPowSpec(:,306,:) = [];
exp_z(306,:) = [];
offs_z(306,:) = [];


%% DELTA
subj=[];
parcel=[];
age=[];
feat=[];
for k1=1:349 %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    age = [age ones(1,200)*age_years(k1)];
    feat = [feat squeeze(corrPowSpec(1,k1,:))'];
end

for rot  = 1:10000
    if ~mod(rot,1000)
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
xlabel('LMEM t-values for permutated StatAvl100 map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(60.66,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_corrdelta_tvalues_thickness.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
p(2,1) =  (10000-max(find(lm_sorted<60.66)))/10000;


%% THETA
subj=[];
parcel=[];
age=[];
feat=[];
for k1=1:349 %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    age = [age ones(1,200)*age_years(k1)];
    feat = [feat squeeze(corrPowSpec(2,k1,:))'];
end

for rot  = 1:45346
    if ~mod(rot,1000)
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
xlabel('LMEM t-values for permutated StatAvl100 map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(48.71,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_corrtheta_tvalues_thickness.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
p(2,2) =  (10000-max(find(lm_sorted<48.71)))/10000;


%% ALPHA
subj=[];
parcel=[];
age=[];
feat=[];
for k1=1:349 %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    age = [age ones(1,200)*age_years(k1)];
    feat = [feat squeeze(corrPowSpec(3,k1,:))'];
end

for rot  = 1:45346
    if ~mod(rot,1000)
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
xlabel('LMEM t-values for permutated StatAvl100 map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(-45.89,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_corralpha_tvalues_thickness.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
p(2,3) =  (10000-max(find(lm_sorted<-45.89)))/10000;


%% BETA
subj=[];
parcel=[];
age=[];
feat=[];
for k1=1:349 %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    age = [age ones(1,200)*age_years(k1)];
    feat = [feat squeeze(corrPowSpec(4,k1,:))'];
end

for rot  = 1:45346
    if ~mod(rot,1000)
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
xlabel('LMEM t-values for permutated StatAvl100 map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(-0.13,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_corrbeta_tvalues_thickness.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
p(2,4) =  (10000-max(find(lm_sorted<-0.13)))/10000;


%% GAMMA
subj=[];
parcel=[];
age=[];
feat=[];
for k1=1:349 %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    age = [age ones(1,200)*age_years(k1)];
    feat = [feat squeeze(corrPowSpec(5,k1,:))'];
end

for rot  = 1:45346
    if ~mod(rot,1000)
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
xlabel('LMEM t-values for permutated StatAvl100 map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(37.62,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_corrgamma_tvalues_thickness.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
p(2,5) =  (10000-max(find(lm_sorted<37.62)))/10000;


%% FOOOF EXPONENT
subj=[];
parcel=[];
age=[];
feat=[];
for k1=1:349 %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    age = [age ones(1,200)*age_years(k1)];
    feat = [feat exp_z(k1,:)];
end

for rot  = 1:10000
    if ~mod(rot,1000)
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
xlabel('LMEM t-values for permutated StatAvl100 map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(-21.99,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_fooofexp_tvalues_thickness.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
p(2,6) =  (10000-max(find(lm_sorted<-21.99)))/10000;


%% FOOOF OOFSET
subj=[];
parcel=[];
age=[];
feat=[];
for k1=1:349 %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    age = [age ones(1,200)*age_years(k1)];
    feat = [feat offs_z(k1,:)'];
end

for rot  = 1:45346
    if ~mod(rot,1000)
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
xlabel('LMEM t-values for permutated StatAvl100 map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];
hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)
xline(-46.70,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100
box off
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_fooofoffs_tvalues_thickness.png'],'Resolution',500)

% which p-value does the t-vual correspond to?
p(2,7) =  (10000-max(find(lm_sorted<-46.70)))/10000;








