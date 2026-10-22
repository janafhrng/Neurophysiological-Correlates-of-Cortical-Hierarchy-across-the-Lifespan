clearvars;
close all
clc

filepath = pwd;

%loading matrix with the data in the form of features x parcels
% matrix is already z-scored 
load([filepath '\DATA\mat_z.mat'])

% loading the feature file with the complete list of features 
features = readtable([filepath '\DATA\features_z.csv']);

% load demographics. columns: participant code,age,sex(1 = male),
% intraccranial volume
demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));

% load the parcellated map
load([filepath '\DATA\hierarchy_parc.csv'])

% extract the statistics of the linear mixed effect model for each feature
subj=[];
parcel=[];
feat=[];
age=[];
hierarchy=[];

for k1=1:size(mat_z,1) %subjects
    subj=[subj ones(1,200)*k1];
    parcel=[parcel 1:200];
    age = [age ones(1,200)*age_years(k1)];
    hierarchy=[hierarchy hierarchy_parc'];
end

for ifeature=1:size(mat_z,3) % loop over features
    if (mod(ifeature,200)==0)
        disp(ifeature);
    end
    feat=[];
    for k1=1:size(mat_z,1) %subjects
        feat=[feat squeeze(mat_z(k1,:,ifeature))];
    end
    varnames={'subj','parcel','feat','age','hierarchy'};
    tbl=table(subj',parcel',feat',age',hierarchy','VariableNames',varnames);
    model = 'feat ~ hierarchy*age + (1|subj) + (1|parcel)';
    try
        lme=fitlme(tbl,model);
        allst(ifeature,:)=lme.Coefficients.tStat;
        pval(ifeature,:)=lme.Coefficients.pValue;
    end
end

% save lme statistics
save([filepath '\DATA\lme_out_age_effect_hierarchy.mat'],'allst','pval', '-v7.3');
writematrix(allst,[filepath '\DATA\lme_out_age_effect_hierarchy.csv']);



%% PLOT FOR SUPPLEMENTARY

load([filepath '\DATA\lme_out_age_effect_hierarchy.mat'])

% plot distribution of t-values main effect hierarchy
index = find(abs(allst(:,3)));
figure('Color','white')
hist(abs(allst(index,3)),50)
h = findobj(gca,'Type','patch');
h.FaceColor = [0 0.5 0.5];
h.FaceAlpha = 0.7;
h.EdgeColor = 'w';
h.EdgeAlpha = 0.1;
set(gca,'xcolor','black','ycolor','black')
xlabel('absolute t statistics','FontSize',16)
ylabel('count','FontSize',16)
%ylim([-0.15,40])
%xlim([69.9 90])
set(gca,'FontSize',14)
box off
hold on;
%line([79.5 79.5], [0 39],'LineStyle','--', 'LineWidth', 2, 'Color',[.7 .7 .7]);
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_main_hierarchy_tValues.png'],'Resolution',500)

% plot distribution of t-values interaction effect hierarchy and age
index = find(abs(allst(:,4)));
figure('Color','white')
hist(abs(allst(index,4)),50)
h = findobj(gca,'Type','patch');
h.FaceColor = [0 0.5 0.5];
h.FaceAlpha = 0.7;
h.EdgeColor = 'w';
h.EdgeAlpha = 0.1;
set(gca,'xcolor','black','ycolor','black')
xlabel('absolute t statistics','FontSize',16)
ylabel('count','FontSize',16)
%ylim([-0.2,90])
%xlim([19.95 55])
set(gca,'FontSize',14)
box off
hold on;
%line([31 31], [0 85],'LineStyle','--', 'LineWidth', 2, 'Color',[.7 .7 .7]);
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_hist_interaction_hierarchy_age_tValues.png'],'Resolution',500)

