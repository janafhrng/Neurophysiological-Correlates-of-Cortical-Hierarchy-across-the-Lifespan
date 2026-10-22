%% REVIEWER 1: if we do multiple linear regression, which features from the whole set predict cortical hierarchy best?
% lasso!!

clearvars;
close all
clc

filepath = pwd;

%loading matrix with the data in the form of features x parcels
% matrix is already z-scored 
load([filepath '\DATA\mat_z.mat'])
% average across subjects
mat_avg = squeeze(mean(mat_z,1));

% loading the feature file with the complete list of features 
features = readtable([filepath '\DATA\features_z.csv']);

% load demographics. columns: participant code,age,sex(1 = male),
% intraccranial volume
demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));
sex = table2array(demographics(:,3)); % sex(1 = female)

% load the parcellated map
load([filepath '\DATA\hierarchy_parc.csv'])

% LASSO REGRESSION IN A NESTED CROSS-VALIDATION LOOP WITH 50 RUNS

for iruns = 1:50
    cv = cvpartition(200,"Kfold",10);

    for ifold = 1:10
        trainIdx(ifold,:) = training(cv, ifold);
        testIdx(ifold,:) = test(cv, ifold);

        %compute lasso
        [B,FitInfo] = lasso(mat_avg(trainIdx(ifold,:),:),hierarchy_parc(trainIdx(ifold,:)),'CV', 10);

        % save paramenters
        idxLambda1SE = FitInfo.Index1SE;
        coef = B(:,idxLambda1SE);
        beta(iruns,ifold,:) = coef;
        coef0 = FitInfo.Intercept(idxLambda1SE);
        interc(iruns,ifold) = coef0;

        % now predict in test data
        tmp_test=squeeze(mat_avg(testIdx(ifold,:),:));
        Y1(cv.test(ifold))= tmp_test * coef + coef0;

        clear coef coef0
    end

    % save prediction performance for each repetition
    acc1(1,iruns)=corr(Y1',hierarchy_parc); %compute correlation
    mae1(1,iruns)=mean(abs(Y1'-hierarchy_parc)); %compute mean absolute error
    mean_obs=mean(hierarchy_parc);
    r2_1(1,iruns)=1-(mean((hierarchy_parc-Y1').^2)/(mean((hierarchy_parc-mean_obs).^2))); %compute R2 score
    Y_pred(:,iruns) = Y1;
end

% average across repetitions
av_acc = mean(acc1, 2);
av_mae = mean(mae1, 2);
av_r2 = mean(r2_1,2);
av_yhat = mean(Y_pred, 2);


% now, look coser into the betas
% average betas across folds, then across repetitions
for iruns = 1:50
 av_beta_folds(iruns,:) = mean(abs(beta(iruns,:,:)), 2);
 av_beta_folds_notabs(iruns,:) = mean(beta(iruns,:,:), 2);
end

av_allbeta = mean(av_beta_folds, 1);
av_allbeta_notabs = mean(av_beta_folds_notabs, 1);


figure('Color','white','Position',[1,1,700,500])
hold on
scatter(hierarchy_parc,av_yhat,90,vik(40,:),'filled','o','MarkerFaceAlpha',.5)
plot(hierarchy_parc,hierarchy_parc,'LineWidth',3,'Color',[0.5 0.5 0.5])
xlabel('Actual Cortical Hierarchy','FontSize',16)
ylabel('Predicted Cortical Hierarchy','FontSize',16)
set(gca,'box','off') 
set(gca,'FontSize',20)
exportgraphics(gcf,[filepath '\PLOTS\SUPPL_lasso_acutal_predicted_hierarchy.png'],'Resolution',500)


% % eport for r??
tmp_idx = find(av_allbeta ~= 0);
idx_mat = mat_avg(:,tmp_idx);
idx_feature = features(tmp_idx,:);
coef = av_allbeta(tmp_idx);
writematrix(idx_mat,[filepath '\DATA\lasso_mat.csv']);
writetable(idx_feature,[filepath '\DATA\lasso_feat.csv'],'WriteRowNames',true);
writematrix(coef,[filepath '\DATA\lasso_coef.csv']);


%% PLOT THE PREDICTED HIERARCHY
[surf_lh, surf_rh] = fetch_template_surface('fsaverage5');
schaefer_200 = fetch_parcellation('fsaverage5', 'schaefer', 200);
load vik.mat

obj = plot_hemispheres(av_yhat, {surf_lh,surf_rh}, ...
    'parcellation',schaefer_200, 'views', 'lm');
set(gcf,'units','points','position',[50,-50,1200,1200])
colormap([0.7 0.7	0.7;vik])
saveas(gcf,['x.png'])
images = imcrop(imread(['x.png']),[170, 30, 1105, 270]);
imwrite(images,[filepath '\PLOTS\SUPPL_predictedHierarchyLasso.png']);

% %% plot the 5 features with the highest weights!!
% [beta_sort,idx_max] = sort(abs(av_allbeta),'descend');
% 
% features(idx_max(1:45),:)
% beta_sort(1:45)
% r_sort = corr(mat_avg(:,idx_max(1:45)),hierarchy_parc)
% 
% results = [features(idx_max(1:45),:), num2cell(beta_sort(1:45))',num2cell(r_sort)];
% 
% %addpath(genpath('C:\Users\fehring\sciebo\Toolboxes'))
% [surf_lh, surf_rh] = load_conte69();
% labeling = load_parcellation('schaefer',200);
% load("vik.mat")
% 
% obj = plot_hemispheres(mat_avg(:,idx_max(1:5)), {surf_lh,surf_rh}, ...
%     'parcellation',labeling.schaefer_200, 'views', 'lm');
% set(gcf,'units','points','position',[50,-50,1200,1200])
% colormap([0.7 0.7	0.7;vik])
% 
% 
% 
% %% are the correlations high between cortical hierarchy and the selected features?
% 
% % spearman
% r = corr(hierarchy_parc,idx_mat,'type','Spearman');
% hist(r)
% 
% % how high is the correaltion with corticla hierarchy for the features with the highest weights?
% r(idx_max5(1:5))
% 
% % which are the features with the highest correlation?
% [r_sort,idx_r] = sort(abs(r),'descend');
% idx_feature(idx_r(1:10),:)
% 
% obj = plot_hemispheres(idx_mat(:,idx_r(1:5)), {surf_lh,surf_rh}, ...
%     'parcellation',labeling.schaefer_200, 'views', 'lm');
% set(gcf,'units','points','position',[50,-50,1200,1200])
% colormap([0.7 0.7	0.7;vik])
% 
% r(idx_r(1:5))
% coef(idx_r(1:5))






