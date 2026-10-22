filepath = pwd;

%% what kind of t-values can be expected?

%% STATAVL100

% use the permutations from python where we create 50.000 permutations
% with alexander bloch 
permutations = readmatrix([filepath '\DATA\spins_StatAvl100_permutations.csv']);

% load demogrpahics in order to catch the age in years
demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));      


% load the parcellated map of cortical hierarchy (sydnor et al., 2021)
load([filepath '\DATA\hierarchy_parc.csv'])


%% null distribution for randomly expactable t-values
% extract the statistics of the linear mixed effect model for each feature
subj=[];
parcel=[];
age=[];
hierarchy=[];
y=[];
for k1=1:size(mat_z,1) %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    age = [age ones(1,200)*age_years(k1)];
    hierarchy=[hierarchy hierarchy_parc'];
end


for rot  = 1:50000
    perm=[];
    for s = 1:350
        perm =[perm permutations(:,rot)'];
    end

    varnames={'subj','parcel','perm','age','hierarchy'};
    tbl=table(subj',parcel',perm',age',hierarchy','VariableNames',varnames);
    model = 'perm ~ hierarchy*age + (1|subj)';

    try
        lme=fitlme(tbl,model);

        %plotResiduals(lme,'fitted')
        lm(rot)=lme.Coefficients.tStat(3);
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*50000;
idxmin1 = 0.005*50000;
idxmax5 = 0.975*50000;
idxmax1 = 0.995*50000;

lm_sorted = sort(lm,"ascend");
lm_sorted(idxmax5) % 62.65
lm_sorted(idxmax1) % 80.2
lm_sorted(idxmin5) % -60.60
lm_sorted(idxmin1) % -74.67

figure('Color','w')
h = histogram(lm,'NumBins',20)
xlabel('LMEM t-values for permutated StatAvl100 map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];

hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)

xline(80.11,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100

box off

exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_perm_statavl_tvalues.png'],'Resolution',500)



%% DELTA
% use the permutations from python where we create 50.000 permutations
% with alexander bloch 
permutations = readmatrix([filepath '\DATA\spins_delta_permutations.csv']);

% extract the statistics of the linear mixed effect model for each feature
subj=[];
parcel=[];
age=[];
hierarchy=[];
y=[];
for k1=1:size(mat_z,1) %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    age = [age ones(1,200)*age_years(k1)];
    hierarchy=[hierarchy hierarchy_parc'];
end


for rot  = 1:50000
    perm=[];
    for s = 1:350
        perm =[perm permutations(:,rot)'];
    end

    varnames={'subj','parcel','perm','age','hierarchy'};
    tbl=table(subj',parcel',perm',age',hierarchy','VariableNames',varnames);
    model = 'perm ~ hierarchy*age + (1|subj)';

    try
        lme=fitlme(tbl,model);

        %plotResiduals(lme,'fitted')
        lm2(rot)=lme.Coefficients.tStat(3);
    end
end

% find indices for 2.5% and 97.5% where the t-valzes should be significant
idxmin5 = 0.025*50000;
idxmin1 = 0.005*50000;
idxmax5 = 0.975*50000;
idxmax1 = 0.995*50000;

lm_sorted = sort(lm2,"ascend");
lm_sorted(idxmax5) % 62.65
lm_sorted(idxmax1) % 80.2
lm_sorted(idxmin5) % -60.60
lm_sorted(idxmin1) % -74.67

figure('Color','w')
h = histogram(lm2,'NumBins',20)
xlabel('LMEM t-values for permutated theta power map')
ylabel('count')
h.FaceColor = [0 0.5 0.5];

hold on 
xline([lm_sorted(idxmax5) lm_sorted(idxmax1) lm_sorted(idxmin5) lm_sorted(idxmin1)],'.','LineWidth',2,'Alpha',0.7)

xline(78.1,'r--','LineWidth',2,'Alpha',0.7) % StatAvl100

box off

exportgraphics(gcf,[filepath '\PLOTS\hist_perm_delta_tvalues.png'],'Resolution',500)




