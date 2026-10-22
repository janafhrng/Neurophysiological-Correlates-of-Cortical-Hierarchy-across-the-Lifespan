%% LIST OF THE 50 TOP FEATURES IN LMEM FOR MAIN CORTICAL HIERARCHY

filepath = pwd;

% loading the feature file with the complete list of features 
features = readtable([filepath '\DATA\features_z.csv']);

% loading lme stats
allst = readmatrix([filepath '\DATA\lme_out_age_effect_hierarchy.csv']);

% ranking the features
% we rank the features according to the highes t-score in order to find
% features that show a really strong gradient along cortical hierarchy
t = allst(:,3);
[~,rank] = sort(abs(t),'descend');
t_ranked = t(rank);
features_ranked = features(rank,:);

data2plot = [features_ranked(1:50,[2,4]),array2table(t_ranked(1:50))];
writetable(data2plot,[filepath '\PLOTS\SUPPL_list_top50_main_corticalHierarchy.csv'])

%% LIST OF THE 50 TOP FEATURES IN LMEM FOR Interaction Effect CORTICAL HIERARCHY * AGE

% ranking the features
% we rank the features according to the highes t-score in order to find
% features that show a really strong gradient along cortical hierarchy *
% age
t = allst(:,4);
[~,rank] = sort(abs(t),'descend');
t_ranked = t(rank);
features_ranked = features(rank,:);

data2plot = [features_ranked(1:50,[2,4]),array2table(t_ranked(1:50))];
writetable(data2plot,[filepath '\PLOTS\SUPPL_list_top50_interaction_corticalHierarchy&Age.csv'])



%% statavl (#448) & hierarchy

%loading matrix with the data in the form of features x parcels
% matrix is already z-scored 
load([filepath '\DATA\mat_z.mat'])

% load demographics. columns: participant code,age,sex(1 = male),
% intraccranial volume
demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));

% load the parcellated map
load([filepath '\DATA\hierarchy_parc.csv'])

% AC40 (#123) & hierarchy
for ifeature=123%:size(mat_z,3) % loop over features
    subj=[];
    parcel=[];
    feat=[];
    age=[];
    hierarchy=[];
    for k1=1:size(mat_z,1) %subjects
        subj=[subj ones(1,200)*k1];
        parcel=[parcel 1:200];
        age=[age ones(1,200)*age_years(k1)];
        feat=[feat squeeze(mat_z(k1,:,ifeature))];
        hierarchy=[hierarchy hierarchy_parc'];
    end
    varnames={'subj','parcel','feat','age','hierarchy'};
    tbl=table(subj',parcel',feat',age',hierarchy','VariableNames',varnames);
    model= 'feat ~ age*hierarchy + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        %plotResiduals(lme,'fitted')
        allst(ifeature,:)=lme.Coefficients.tStat;
        pval(ifeature,:)=lme.Coefficients.pValue;
    end
end

tab1 = table2cell(dataset2table(lme.Coefficients));
xlswrite([filepath '\PLOTS\SUPPL_LMEMresults_AC40_hierarchy'],tab1)


%% statavl (#448) & thickness
mat_z(306,:,:) = [];
% loading cortical thickness
load([filepath '\DATA\thickness_smoothed_12_schaefer200_reordered.mat']);
thickness(306,:) = [];
% load age data
demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));   
age_years(306) = [];


for ifeature=123%1:size(mat_z,3) % loop over features
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
    model= 'feat ~ age*thick + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        %plotResiduals(lme,'fitted')
        allst(ifeature,:)=lme.Coefficients.tStat;
        pval(ifeature,:)=lme.Coefficients.pValue;
    end
end

tab1 = table2cell(dataset2table(lme.Coefficients));
xlswrite([filepath '\PLOTS\SUPPL_LMEMresults_AC40_thickness'],tab1)




%% delta and hierarchy
% loading the the power of different freq bands over all parcels
load([filepath '\DATA\03_specData_freqbands_norm.mat'])
freqbands{1} = zscore(delta_norm,[],2); 
load([filepath '\DATA\hierarchy_parc.csv'])
demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));     

for ifreq = 1%:5 % looping through the canonical frequency bands
    subj=[];
    parcel=[];
    feat=[];
    age=[];
    hierarchy=[];
    for k1=1:350 %subjects
        subj=[subj ones(1,200)*k1];
        parcel=[parcel 1:200];
        age=[age ones(1,200)*age_years(k1)];
        feat=[feat freqbands{1}(k1,:)];
        hierarchy=[hierarchy hierarchy_parc'];
    end
    varnames={'subj','parcel','feat','age','hierarchy'};
    tbl=table(subj',parcel',feat',age',hierarchy','VariableNames',varnames);
    model= 'feat ~ age*hierarchy + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        %plotResiduals(lme,'fitted')
        allst(ifeature,:)=lme.Coefficients.tStat;
        pval(ifeature,:)=lme.Coefficients.pValue;
    end
end

tab1 = table2cell(dataset2table(lme.Coefficients));
xlswrite([filepath '\PLOTS\SUPPL_LMEMresults_delta_hierarchy'],tab1)

%% delta and thickness
% loading the the power of different freq bands over all parcels
freqbands{1}(306,:) = [];
% loading cortical thickness
load([filepath '\DATA\thickness_smoothed_12_schaefer200_reordered.mat']);
thickness(306,:) = [];
age_years(306,:) = [];       

for ifreq = 1%:5 % looping through the canonical frequency bands
    subj=[];
    parcel=[];
    feat=[];
    age=[];
    thick=[];
    for k1=1:349 %subjects
        subj=[subj ones(1,200)*k1];
        parcel=[parcel 1:200];
        age=[age ones(1,200)*age_years(k1)];
        feat=[feat freqbands{1}(k1,:)];
        thick=[thick thickness(k1,:)];
    end
    varnames={'subj','parcel','feat','age','thick'};
    tbl=table(subj',parcel',feat',age',thick','VariableNames',varnames);
    model= 'feat ~ age*thick + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        %plotResiduals(lme,'fitted')
        allst(ifeature,:)=lme.Coefficients.tStat;
        pval(ifeature,:)=lme.Coefficients.pValue;
    end
end

tab1 = table2cell(dataset2table(lme.Coefficients));
xlswrite([filepath '\PLOTS\SUPPL_LMEMresults_delta_thickness'],tab1)


%% ALPHA and HIERARCHY and AGE
% loading the the power of different freq bands over all parcels
load([filepath '\DATA\03_specData_freqbands_norm.mat'])
freqbands{1} = zscore(alpha_norm,[],2); 
load([filepath '\DATA\hierarchy_parc.csv'])
demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));     

for ifreq = 1%:5 % looping through the canonical frequency bands
    subj=[];
    parcel=[];
    feat=[];
    age=[];
    hierarchy=[];
    for k1=1:350 %subjects
        subj=[subj ones(1,200)*k1];
        parcel=[parcel 1:200];
        age=[age ones(1,200)*age_years(k1)];
        feat=[feat freqbands{1}(k1,:)];
        hierarchy=[hierarchy hierarchy_parc'];
    end
    varnames={'subj','parcel','feat','age','hierarchy'};
    tbl=table(subj',parcel',feat',age',hierarchy','VariableNames',varnames);
    model= 'feat ~ age*hierarchy + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        %plotResiduals(lme,'fitted')
        allst(ifeature,:)=lme.Coefficients.tStat;
        pval(ifeature,:)=lme.Coefficients.pValue;
    end
end

tab1 = table2cell(dataset2table(lme.Coefficients));
xlswrite([filepath '\PLOTS\SUPPL_LMEMresults_alpha_hierarchy'],tab1)

%% ALPHA and THICKNESS and AGE
% loading the the power of different freq bands over all parcels
freqbands{1}(306,:) = [];
% loading cortical thickness
load([filepath '\DATA\thickness_smoothed_12_schaefer200_reordered.mat']);
thickness(306,:) = [];
age_years(306,:) = [];       

for ifreq = 1%:5 % looping through the canonical frequency bands
    subj=[];
    parcel=[];
    feat=[];
    age=[];
    thick=[];
    for k1=1:349 %subjects
        subj=[subj ones(1,200)*k1];
        parcel=[parcel 1:200];
        age=[age ones(1,200)*age_years(k1)];
        feat=[feat freqbands{1}(k1,:)];
        thick=[thick thickness(k1,:)];
    end
    varnames={'subj','parcel','feat','age','thick'};
    tbl=table(subj',parcel',feat',age',thick','VariableNames',varnames);
    model= 'feat ~ age*thick + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        %plotResiduals(lme,'fitted')
        allst(ifeature,:)=lme.Coefficients.tStat;
        pval(ifeature,:)=lme.Coefficients.pValue;
    end
end

tab1 = table2cell(dataset2table(lme.Coefficients));
xlswrite([filepath '\PLOTS\SUPPL_LMEMresults_alpha_thickness'],tab1)

%% GAMMA and HIERARCHY and AGE
% loading the the power of different freq bands over all parcels
load([filepath '\DATA\03_specData_freqbands_norm.mat'])
freqbands{1} = zscore(gamma_norm,[],2); 
load([filepath '\DATA\hierarchy_parc.csv'])
demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));     

for ifreq = 1%:5 % looping through the canonical frequency bands
    subj=[];
    parcel=[];
    feat=[];
    age=[];
    hierarchy=[];
    for k1=1:350 %subjects
        subj=[subj ones(1,200)*k1];
        parcel=[parcel 1:200];
        age=[age ones(1,200)*age_years(k1)];
        feat=[feat freqbands{1}(k1,:)];
        hierarchy=[hierarchy hierarchy_parc'];
    end
    varnames={'subj','parcel','feat','age','hierarchy'};
    tbl=table(subj',parcel',feat',age',hierarchy','VariableNames',varnames);
    model= 'feat ~ age*hierarchy + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        %plotResiduals(lme,'fitted')
        allst(ifeature,:)=lme.Coefficients.tStat;
        pval(ifeature,:)=lme.Coefficients.pValue;
    end
end

tab1 = table2cell(dataset2table(lme.Coefficients));
xlswrite([filepath '\PLOTS\SUPPL_LMEMresults_gamma_hierarchy'],tab1)

%% GAMME and THICKNESS and AGE
% loading the the power of different freq bands over all parcels
freqbands{1}(306,:) = [];
% loading cortical thickness
load([filepath '\DATA\thickness_smoothed_12_schaefer200_reordered.mat']);
thickness(306,:) = [];
age_years(306,:) = [];       

for ifreq = 1%:5 % looping through the canonical frequency bands
    subj=[];
    parcel=[];
    feat=[];
    age=[];
    thick=[];
    for k1=1:349 %subjects
        subj=[subj ones(1,200)*k1];
        parcel=[parcel 1:200];
        age=[age ones(1,200)*age_years(k1)];
        feat=[feat freqbands{1}(k1,:)];
        thick=[thick thickness(k1,:)];
    end
    varnames={'subj','parcel','feat','age','thick'};
    tbl=table(subj',parcel',feat',age',thick','VariableNames',varnames);
    model= 'feat ~ age*thick + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        %plotResiduals(lme,'fitted')
        allst(ifeature,:)=lme.Coefficients.tStat;
        pval(ifeature,:)=lme.Coefficients.pValue;
    end
end

tab1 = table2cell(dataset2table(lme.Coefficients));
xlswrite([filepath '\PLOTS\SUPPL_LMEMresults_gamma_thickness'],tab1)


%% 1/f and HIERARCHY and AGE
load([filepath '\DATA\correctedFooofParam_fooofWknee_z.mat']);
load([filepath '\DATA\hierarchy_parc.csv'])
demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));     

subj=[];
parcel=[];
feat=[];
age=[];
hierarchy=[];
for k1=1:350 %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    age=[age ones(1,200)*age_years(k1)];
    feat=[feat exp_z(k1,:)];
    hierarchy=[hierarchy hierarchy_parc'];
end
varnames={'subj','parcel','feat','age','hierarchy'};
tbl=table(subj',parcel',feat',age',hierarchy','VariableNames',varnames);
model= 'feat ~ age*hierarchy + (1|subj) + (1|parcel)';
try
    lme=fitlme(tbl,model);
    %plotResiduals(lme,'fitted')
    allst(ifeature,:)=lme.Coefficients.tStat;
    pval(ifeature,:)=lme.Coefficients.pValue;
end

tab1 = table2cell(dataset2table(lme.Coefficients));
xlswrite([filepath '\PLOTS\SUPPL_LMEMresults_aperidocExp_hierarchy'],tab1)

%% 1/f and THICKNESS and AGE
exp_z(306,:) = [];
% loading cortical thickness
load([filepath '\DATA\thickness_smoothed_12_schaefer200_reordered.mat']);
thickness(306,:) = [];
age_years(306,:) = [];       

subj=[];
parcel=[];
feat=[];
age=[];
thick=[];
for k1=1:349 %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    age=[age ones(1,200)*age_years(k1)];
    feat=[feat exp_z(k1,:)];
    thick=[thick thickness(k1,:)];
end
varnames={'subj','parcel','feat','age','thick'};
tbl=table(subj',parcel',feat',age',thick','VariableNames',varnames);
model= 'feat ~ age*thick + (1|subj) + (1|parcel)';
try
    lme=fitlme(tbl,model);
    %plotResiduals(lme,'fitted')
    allst(ifeature,:)=lme.Coefficients.tStat;
    pval(ifeature,:)=lme.Coefficients.pValue;
end

tab1 = table2cell(dataset2table(lme.Coefficients));
xlswrite([filepath '\PLOTS\SUPPL_LMEMresults_aperiodicExp_thickness'],tab1)



%% AC40 & hierarchy & AGE

%loading matrix with the data in the form of features x parcels
% matrix is already z-scored 
load([filepath '\DATA\mat_z.mat'])

% load demographics. columns: participant code,age,sex(1 = male),
% intraccranial volume
demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));

% load the parcellated map
load([filepath '\DATA\hierarchy_parc.csv'])

% AC40 (#123) & hierarchy
for ifeature=123%:size(mat_z,3) % loop over features
    subj=[];
    parcel=[];
    feat=[];
    age=[];
    hierarchy=[];
    for k1=1:size(mat_z,1) %subjects
        subj=[subj ones(1,200)*k1];
        parcel=[parcel 1:200];
        age=[age ones(1,200)*age_years(k1)];
        feat=[feat squeeze(mat_z(k1,:,ifeature))];
        hierarchy=[hierarchy hierarchy_parc'];
    end
    varnames={'subj','parcel','feat','age','hierarchy'};
    tbl=table(subj',parcel',feat',age',hierarchy','VariableNames',varnames);
    model= 'feat ~ age*hierarchy + (1|subj)';
    try
        lme=fitlme(tbl,model);
        %plotResiduals(lme,'fitted')
        allst(ifeature,:)=lme.Coefficients.tStat;
        pval(ifeature,:)=lme.Coefficients.pValue;
    end
end

tab1 = table2cell(dataset2table(lme.Coefficients));
xlswrite([filepath '\PLOTS\SUPPL_LMEMresults_AC40_hierarchy'],tab1)


%% Ac40 & thickness & age
mat_z(306,:,:) = [];
% loading cortical thickness
load([filepath '\DATA\Input\thickness_smoothed_12_schaefer200_reordered.mat']);
thickness(306,:) = [];
% load age data
demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));   
age_years(306) = [];


for ifeature=123%1:size(mat_z,3) % loop over features
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
    model= 'feat ~ age*thick + (1|subj)';
    try
        lme=fitlme(tbl,model);
        %plotResiduals(lme,'fitted')
        allst(ifeature,:)=lme.Coefficients.tStat;
        pval(ifeature,:)=lme.Coefficients.pValue;
    end
end

tab1 = table2cell(dataset2table(lme.Coefficients));
xlswrite([filepath '\PLOTS\SUPPL_LMEMresults_AC40_thickness'],tab1)


%% AMI40 & hierarchy & AGE

% save ami
load([filepath '\DATA\ami_1_40.mat'])

% load demographics. columns: participant code,age,sex(1 = male),
% intraccranial volume
demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));

% load the parcellated map
load([filepath '\DATA\hierarchy_parc.csv'])

% AMI40 (#40) & hierarchy
for ifeature=40
    subj=[];
    parcel=[];
    feat=[];
    age=[];
    hierarchy=[];
    for k1=1:350 %subjects
        subj=[subj ones(1,200)*k1];
        parcel=[parcel 1:200];
        age=[age ones(1,200)*age_years(k1)];
        feat=[feat squeeze(ami_out(k1,:,ifeature))];
        hierarchy=[hierarchy hierarchy_parc'];
    end
    varnames={'subj','parcel','feat','age','hierarchy'};
    tbl=table(subj',parcel',feat',age',hierarchy','VariableNames',varnames);
    model= 'feat ~ age*hierarchy + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        %plotResiduals(lme,'fitted')
        allst(ifeature,:)=lme.Coefficients.tStat;
        pval(ifeature,:)=lme.Coefficients.pValue;
    end
end

tab1 = table2cell(dataset2table(lme.Coefficients));
xlswrite([filepath '\PLOTS\SUPPL_LMEMresults_AMI40_hierarchy'],tab1)


%% AMI40 & thickness & age
ami_out(306,:,:) = [];
% loading cortical thickness
load([filepath '\DATA\thickness_smoothed_12_schaefer200_reordered.mat']);
thickness(306,:) = [];
% load age data
demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));   
age_years(306) = [];


for ifeature=40%1:size(mat_z,3) % loop over features
    subj=[];
    parcel=[];
    thick=[];
    feat=[];
    age=[];
    for k1=1:349 %subjects
        subj=[subj ones(1,200)*k1];
        parcel=[parcel 1:200];
        thick=[thick thickness(k1,:)];
        feat=[feat squeeze(ami_out(k1,:,ifeature))];
        age = [age ones(1,200)*age_years(k1)];
    end
    varnames={'subj','parcel','feat','thick','age'};
    tbl=table(subj',parcel',feat',thick',age','VariableNames',varnames);
    model= 'feat ~ age*thick + (1|subj)+(1|parcel)';
    try
        lme=fitlme(tbl,model);
        %plotResiduals(lme,'fitted')
        allst(ifeature,:)=lme.Coefficients.tStat;
        pval(ifeature,:)=lme.Coefficients.pValue;
    end
end

tab1 = table2cell(dataset2table(lme.Coefficients));
xlswrite([filepath '\PLOTS\SUPPL_LMEMresults_AMI40_thickness'],tab1)
