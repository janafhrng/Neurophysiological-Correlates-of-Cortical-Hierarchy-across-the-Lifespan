filepath = 'G:\JF\HCTSA feature gradients\CamCan\HierarchyCorrelates';

%% Alpha POWER ÜBER GRADIOMETER UND MAGNETOMETER
% load u.a. fieldtrip toolbox

% load time series data files
% data: for each participant: parcels(214) x trials(30) x samples(3000)
PATHIN = 'C:\Users\fehringj\Sciebo\TBraic\JF\all_subjects_lcmv_grads\';
% only look for files which are from condition 1 (closed eyes condition)
files = dir([PATHIN,'*.mat']);
files = {files.name};

load('G:\JF\HCTSA feature gradients\CamCan\HierarchyCorrelates\DATA\mapping_idx_camcan2schaefer.mat')


%% ALPHA POWER!!! ========================================================
for i = 1:size(files,2)
    i
    % load the specific data file
    load([PATHIN files{i}]);

    % create fieldtrip data struct
    data=[];
    data.fsample=300;
    for k=1:214
        data.label{k}=num2str(k);
    end
    for k=1:30
        data.trial{k} = squeeze(ts(:,k,:));
        data.time{k} = 0:1/300:9.999;
    end
    data.dimord='chan_time';
    
    % setup configuration
    cfg = [];
    cfg.output = 'pow';
    cfg.method = 'mtmfft'; % multitaper frequency transformation
    cfg.taper = 'dpss';
    cfg.foi = [0:100]; % complete band
    cfg.tapsmofrq = 1;
    cfg.pad = 15;
    %cfg.keeptrials = 'yes';

    % do freq analysis
    tmp_specData = ft_freqanalysis(cfg, data);
    % average over the 6 frequencies
    specData_all(i,:,:) = tmp_specData.powspctrm;
end

specData = specData_all(:,tmp_idx,:);

save('G:\JF\HCTSA feature gradients\CamCan\HierarchyCorrelates\DATA\powerspec_gradiometer.mat','specData')

all = mean(specData,3,'omitnan');
alpha =  mean(specData(:,:,9:14),3,'omitnan');
alpha_norm = alpha./all; 
alpha_norm_z = zscore(alpha_norm,[],2); % zscore 
alpha_norm_mean = mean(alpha_norm_z,1); % mean

% correlate with magnetometer data
clear delta_norm
load('G:\JF\HCTSA feature gradients\CamCan\HierarchyCorrelates\DATA\03_specData_freqbands_norm.mat')
alpha_mag = zscore(alpha_norm,[],2);
alpha_mag_z = zscore(alpha_mag,[],2); % zscore 
alpha_mag_mean = mean(alpha_mag_z,1);


% PLOT -----------------------------------------------------------
parcellation = fetch_parcellation('fslr32k', 'schaefer', 200);
[surf_lh, surf_rh] = fetch_template_surface('fslr32k', 'layer', 'inflated');
load("vik.mat")

data2plot2 = [alpha_norm_mean',alpha_mag_mean'];

vertex_data = nan(size(parcellation,1),2);
for i_sensor = 1:2
for i = 1:200
    vertex_data(parcellation == i,i_sensor) = data2plot2(i,i_sensor);
end
end

obj = plot_hemispheres(vertex_data, {surf_lh, surf_rh});
set(gcf,'units','points','position',[65,-65,1215,1215])
colormap([0.7 0.7 0.7; vik])
%colorlimits(obj,[min(vertex_data,[],'all'), max(vertex_data,[],'all')])
saveas(gcf,'x.png')
I1 = imread('x.png');
saveas(gcf,'x.png');
I1 = imread('x.png');
T1 = imcrop(I1, [170, 50, 1300, 220]); % only left hemisphere
T2 = imcrop(I1, [170, 350, 1300, 220]);
Imerged = imtile({T1;T2},  'BorderSize', 4, 'BackgroundColor', 'white','GridSize', [3,1]);
imwrite(Imerged,'G:\JF\HCTSA feature gradients\CamCan\HierarchyCorrelates\\PLOTS\SUPPL_alpha_power_MAGvsGRAD.png');


% CPRRELATE
corr(alpha_norm_mean', alpha_mag_mean')

  


%% AUTOCORRELATION!! ========================================================

for  isubj = 1:size(files,2)
    isubj
    load([PATHIN files{isubj}]);

    % create fieldtrip data struct
    data=[];
    data.fsample=300;
    for k=1:214
        data.label{k}=num2str(k);
    end
    for k=1:30
        data.trial{k}=squeeze(ts(:,k,:));
        data.time{k}=[0:1/300:9.999];
    end
    data.dimord='chan_time';

    % calculate AC40 across gradiometers
    for itrial = 1:30
        for iparcel = 1:214
            ac40(isubj,itrial,iparcel) = CO_AutoCorr(data.trial{itrial}(iparcel,:),40,'Fourier');
        end
    end
    % average across trials
    statavl100_mean(isubj,:) = squeeze(mean(ac40(isubj,:,:),2))';
    % take only 200 parcels
    load('G:\JF\HCTSA feature gradients\CamCan\HierarchyCorrelates\DATA\mapping_idx_camcan2schaefer.mat')
    AC40_tmp(isubj,:) = statavl100_mean(isubj,tmp_idx);
end

% zscore
AC40_z = zscore(AC40_tmp,[],2); % zscore
%average across all participants
AC40_mean = mean(AC40_z,1);

% compare with VarM magnetometer also only for the first 20 participants
%load([filepath, '\big DATA files\mat_z.mat'])
load([filepath, '\DATA\features_selected.mat'])
% AC40 = #123

AC40_mag = squeeze(mean(mat_z(:,:,123),1));

% PLOT -----------------------------------------------------------------
parcellation = fetch_parcellation('fslr32k', 'schaefer', 200);
[surf_lh, surf_rh] = fetch_template_surface('fslr32k', 'layer', 'inflated');
load("vik.mat")

data2plot2 = [AC40_mean',AC40_mag'];

vertex_data = nan(size(parcellation,1),2);
for i_sensor = 1:2
for i = 1:200
    vertex_data(parcellation == i,i_sensor) = data2plot2(i,i_sensor);
end
end

obj = plot_hemispheres(vertex_data, {surf_lh, surf_rh});
set(gcf,'units','points','position',[65,-65,1215,1215])
colormap([0.7 0.7 0.7; vik])
%colorlimits(obj,[min(vertex_data,[],'all'), max(vertex_data,[],'all')])
saveas(gcf,'x.png')
I1 = imread('x.png');
saveas(gcf,'x.png')
I1 = imread('x.png');
T1 = imcrop(I1, [170, 50, 1300, 220]); % only left hemisphere
T2 = imcrop(I1, [170, 350, 1300, 220]);
Imerged = imtile({T1;T2},  'BorderSize', 4, 'BackgroundColor', 'white','GridSize', [3,1]);
imwrite(Imerged,[filepath '\PLOTS\SUPPL_AC40_MAGvsGRAD.png']);

corr(AC40_mean',AC40_mag')




%% ============================================================================================================================================= 
% re compute LMEM with both features, to see whether these features show a similar dependency on cortical hierarchy/thickness as with magnetometers
load([filepath '\DATA\hierarchy_parc.csv'])

demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));     

features(:,:,1) = alpha_norm_z;
features(:,:,2) = AC40_z;

for ifeat = 1:2 % looping through the two features we are investigatign
    subj=[];
    parcel=[];
    feat=[];
    hierarchy=[];
    age=[];
    for k1=1:350 %subjects
        subj=[subj ones(1,200)*k1];
        feat=[feat features(k1,:,ifeat)];
        hierarchy=[hierarchy hierarchy_parc'];
        age = [age ones(1,200)*age_years(k1)];
        parcel=[parcel 1:200];
    end
    varnames={'subj','feat','hierarchy','age','parcel'};
    tbl=table(subj',feat',hierarchy',age',parcel','VariableNames',varnames);
    model= 'feat ~ age*hierarchy + (1|subj) + (1|parcel)';

    try
        lme = fitlme(tbl,model);
        allst_freqbands(ifeat,:)=lme.Coefficients.tStat;
        pval_freqbands(ifeat,:)=lme.Coefficients.pValue;
    end
end


%% thickness
% loading cortical thickness
load([filepath '\DATA\thickness_smoothed_12_schaefer200_reordered.mat']);
thickness(306,:) = [];
age_years(306,:) = [];   
features(306,:,:) = [];


for ifeat = 1:2 % looping through the two features we are investigatign
    subj=[];
    parcel=[];
    feat=[];
    thick=[];
    age=[];
    for k1=1:349 %subjects
        subj=[subj ones(1,200)*k1];
        feat=[feat features(k1,:,ifeat)];
        thick=[thick thickness(k1,:)];
        age = [age ones(1,200)*age_years(k1)];
        parcel=[parcel 1:200];
    end
    varnames={'subj','feat','thick','age','parcel'};
    tbl=table(subj',feat',thick',age',parcel','VariableNames',varnames);
    model= 'feat ~ age*thick + (1|subj) + (1|parcel)';

    try
        lme = fitlme(tbl,model);
        allst_freqbands(ifeat,:)=lme.Coefficients.tStat;
        pval_freqbands(ifeat,:)=lme.Coefficients.pValue;
    end
end




















