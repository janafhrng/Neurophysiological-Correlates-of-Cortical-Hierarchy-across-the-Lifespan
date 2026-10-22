%% Prepare features for spin with cortical hierarchy

filepath = pwd;
 
% POWER SPECTRUM
% load power spectrum
load([filepath '\DATA\01_specData_all.mat']);

% only take the 200 cortical parcels and disregard the 14 other parcels
load([filepath '\DATA\mapping_idx_camcan2schaefer.mat'])
specData_all = specData_all(:,tmp_idx,:);

% mean of all freqbands
all = mean(specData_all,3,'omitnan');
delta =  mean(specData_all(:,:,2:4),3,'omitnan');
theta =  mean(specData_all(:,:,5:8),3,'omitnan');
alpha =  mean(specData_all(:,:,9:14),3,'omitnan');
beta =  mean(specData_all(:,:,15:31),3,'omitnan');
gamma =  mean(specData_all(:,:,32:101),3,'omitnan');

% normalize over all freqbands
delta_norm = delta./all; 
theta_norm = theta./all; 
alpha_norm = alpha./all; 
beta_norm = beta./all; 
gamma_norm = gamma./all; 

% z-score
delta_z = zscore(delta_norm,[],2);
theta_z = zscore(theta_norm,[],2);
alpha_z = zscore(alpha_norm,[],2);
beta_z = zscore(beta_norm,[],2);
gamma_z = zscore(gamma_norm,[],2);

% take the mean of the freqbands
delta_z_mean = mean(delta_z,1);
theta_z_mean = mean(theta_z,1);
alpha_z_mean = mean(alpha_z,1);
beta_z_mean = mean(beta_z,1);
gamma_z_mean = mean(gamma_z,1);

% FOOOF AND COG
load([filepath '\DATA\foof_variables.mat'])
fooof_z_mean = mean(zscore(fooof_exp_2_60Hz,[],2),1);
% we look at the fooof exponent from 2-60 Hz, the other one is from 1-60 Hz

load([filepath '\DATA\peak_freq_cog.mat'])
peakfreq_z_mean = mean(zscore(peak_freq_cog,[],2),1);

% zscore StatAvl100
%loading matrix with the data in the form of participant x parcels x features
load([filepath '\DATA\all_subj_feat_reordered_schafer200.mat']);
statavl = squeeze(mat(:,:,448)); % statavl100 = feature number 448
statavl_z_mean = mean(zscore(statavl,[],2),1);

% save all features for spin test
feat4spin_hierarchy = [statavl_z_mean; delta_z_mean; theta_z_mean; alpha_z_mean; beta_z_mean; gamma_z_mean; fooof_z_mean; peakfreq_z_mean];
writematrix(feat4spin_hierarchy, [filepath '\DATA\feat4spin_hierarchy.csv'])



%% Prepare features for spin with cortical thickness
% here we need to take one participant out, where the CT data is a large
% outlier (#306)

% z-score
delta_z(306,:) = [];
theta_z(306,:) = [];
alpha_z(306,:) = [];
beta_z(306,:) = [];
gamma_z(306,:) = [];

% take the mean of the freqbands
delta_z_mean = mean(delta_z,1);
theta_z_mean = mean(theta_z,1);
alpha_z_mean = mean(alpha_z,1);
beta_z_mean = mean(beta_z,1);
gamma_z_mean = mean(gamma_z,1);

% FOOOF AND COG
load([filepath '\DATA\foof_variables.mat'])
fooof_exp_2_60Hz(306,:) = []; 
fooof_z_mean = mean(zscore(fooof_exp_2_60Hz,[],2),1);
% we look at the fooof exponent from 2-60 Hz, the other one is from 1-60 Hz

load([filepath '\DATA\peak_freq_cog.mat'])
peak_freq_cog(306,:) = [];
peakfreq_z_mean = mean(zscore(peak_freq_cog,[],2),1);

% zscore StatAvl100
statavl = squeeze(mat(:,:,448)); % statavl100 = feature number 448
statavl(306,:) = [];
statavl_z_mean = mean(zscore(statavl,[],2),1);

% save all features for spin test
feat4spin_hierarchy = [statavl_z_mean; delta_z_mean; theta_z_mean; alpha_z_mean; beta_z_mean; gamma_z_mean; fooof_z_mean; peakfreq_z_mean];
writematrix(feat4spin_hierarchy, [filepath '\DATA\feat4spin_thickness.csv'])


%% PREPARE CORTICAL THICKNESS
% loading cortical thickness
load([filepath '\DATA\thickness_smoothed_12_schaefer200_reordered.mat']);

% take out participant number 306 --> does not have a high correlation with
% the rest of the cortical thickness data
thickness(306,:) = [];
thickness_mean = mean(thickness,1);

writematrix(thickness_mean,[filepath '\DATA\thickness_mean4spin.csv'])











