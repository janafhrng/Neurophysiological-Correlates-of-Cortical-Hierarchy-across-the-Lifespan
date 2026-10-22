%% BONFERRONI CORRECTION

% number of features tested: hctsa: 5961, 7 'traditional' features
numFeat = 5961 + 7;
pBonnf = 0.05/(numFeat*2);

%% HIERARCHY
% load all p-vals from the tested hctsa feautres
filepath = pwd;
load([filepath '\DATA\lme_out_age_effect_hierarchy.mat']);

% test for main effect of hierarchy 
% sort t-values for all 50 highest t-values

% --> sort p-value accordingly!
% --> set 0 pvale to min
% sets
[~,idx] =  sort(abs(allst(:,3)),'descend');
allst_sort = allst(idx,3);
pval_sort = pval(idx,3);
isSign = pval_sort < pBonnf;











