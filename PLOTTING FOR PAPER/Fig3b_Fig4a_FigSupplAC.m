%% PLOTTING AUTOCORRELATION OVER DIFFERENT AGES
filepath = pwd;

load('vik.mat')
load('lipari10.mat')
load('lajolla10.mat')
load("batlow10.mat")

% finding participants in different age ranges
demographics = readtable([filepath '\DATA\demo_all_subjects.csv']);
age_years = table2array(demographics(:,2));
index_age_group = [1:50;151:200;301:350];

% finding parcels which are low, medium or high in cortical hierarchy
load([filepath '\DATA\hierarchy_parc.csv'])
[~,hierarchy_sorted] = sort(hierarchy_parc,'ascend');
idx_h_min = hierarchy_sorted(1:10);
idx_h_med = hierarchy_sorted(96:105);
idx_h_max = hierarchy_sorted(191:200);
hierarchy_labels = {'idx_h_min','idx_h_med','idx_h_max'};

% load feature of interest (here AC1-40)
%loading matrix with the data in the form of participant x parcels x features
load([filepath '\DATA\all_subj_feat_reordered_schafer200.mat']);

% find the features AC1 - AC40
idx = 84:123;
% load the tstatistics for AC1:40
allst = readmatrix([filepath '\DATA\lme_out_age_effect_hierarchy.csv']);
tAC = allst(idx,:);


%% Figure 4a 
colors = [lajolla10(4,:),.8;lajolla10(6,:),.8;lajolla10(8,:),.8];
colorbw = [.9 .9 .9; .6 .6 .6; .3 .3 .3];

for i = 1:5; ticks{i} = num2str(round(i*7.5/300*1000));end

figure('Position', [10 10 350 1000])
for age = 1:3
    ts = mat(index_age_group(age,:),:,:);
    ts1 = squeeze(mean(ts,1));

    % plot the autocorrelation
    subplot(4,1,age)
   
    hold('on');
    for parcel = 1:3
        ts2 = mean(ts1(eval(hierarchy_labels{parcel}),:),1);
        plot(ts2(idx),'-','Color',colors(parcel,:),'LineWidth',3);
    end

    ylim([-0.5, 1.3])
    box off
    xticks([0 7.5 15 22.5 30 37.5])
    xticklabels([0 ticks])
    set(gca,'FontSize',14)
    set(gca, 'Color', 'None');

end

subplot(4,1,4)
for cond = 1:3
    hold('on');
    plot(tAC(:,1+cond),'-','Color',colorbw(cond,:),'LineWidth',2);
    hold('off')
    ylim([-60 50])
    xlabel('lag [ms]')
    xticks([0 7.5 15 22.5 30 37.5])
    xticklabels([0 ticks])
    set(gca,'FontSize',14)
    box off
    set(gca, 'Color', 'None');

end  

exportgraphics(gcf,[filepath '\PLOTS\Fig4a_interarction_AC_age.png'],'Resolution',500)


%% Figure 4 plotting the parcels for legend
% welches sind überhaupt die parcels, die ich mir anschaue? 


colors = [lajolla10(4,:);lajolla10(6,:);lajolla10(8,:)];

data2plot = zeros(1,200);
for parcel = 1:3
    data2plot(eval(hierarchy_labels{parcel})) = parcel;
end

%plot
obj = plot_hemispheres(data2plot', {surf_lh,surf_rh},'parcellation',labeling.schaefer_200);
set(gcf,'units','points','position',[50,-50,1200,1200])
colormap([0.7 0.7 0.7; colors])
saveas(gcf,['x.png'])
images = imcrop(imread(['x.png']),[170, 30, 1105, 270]);
imwrite(images,[filepath '\PLOTS\Fig4_topo_parcels_legend.png']);


%% Figure 3b AC lag 40 across age topography
data2plot = squeeze(mat(:,:,123));

 % split according to age group 
 mat1 = mean(data2plot(table2array(demographics(:,2)) <= 28,:),1);
 mat2 = mean(data2plot(table2array(demographics(:,2)) > 48 & table2array(demographics(:,2)) <= 58,:),1);
 mat3 = mean(data2plot(table2array(demographics(:,2)) > 79,:),1);

 for i = 1:3
    data2plot2(i,:) = eval(['mat', num2str(i)]);
 end

vertex_data = nan(size(parcellation,1),3);
for i_age = 1:3
for i = 1:200
    vertex_data(parcellation == i,i_age) = data2plot2(i_age,i);
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
T1 = imcrop(I1, [170, 50, 565, 220]); % only left hemisphere
T2 = imcrop(I1, [170, 350, 565, 220]);
T3 = imcrop(I1, [170, 650, 565, 220]);
Imerged = imtile({T1;T2;T3},  'BorderSize', 4, 'BackgroundColor', 'white','GridSize', [3,1]);
imwrite(Imerged,[filepath '\PLOTS\Fig3b_AC_topo_across_age.png']);


% plot avergae across all ages
% plot topo at AC40 = AC133ms
data2plot = mean(squeeze(mat(:,:,123)),1);

vertex_data = nan(size(parcellation));
for i = 1:200
    vertex_data(parcellation == i) = data2plot(i);
end

obj = plot_hemispheres(vertex_data, {surf_lh, surf_rh});
set(gcf,'units','points','position',[65,-65,1215,1215])
colormap([0.7 0.7 0.7;vik])
saveas(gcf,'x.png')
I1 = imread('x.png');
T1 = imcrop(I1, [180, 60, 555, 210]); % only left hemisphere
Imerged = imtile({T1},  'BorderSize', 4, 'BackgroundColor', 'white','GridSize', [1,1]);
imwrite(Imerged,[filepath '\PLOTS\Fig3b_AC40_topo.png']);


%% Figure 3b -- scatter and lifit acorss age
% plotting it so the interaction effect of age and cortical hierarchy are
% clear

% only take the 50 youngest and oldest participants and plot their AC40
% value over the the different parcels which have diferent hierarchy levels

colors = {repmat(batlow10(5,:),200,1),repmat(batlow10(7,:),200,1),repmat(batlow10(9,:),200,1)};

for age = [1:3]
    %our feature AC40: 123
    ts = mat(index_age_group(age,:),:,123);
    ts_mean = squeeze(mean(ts,1));

    % plot a linear fit into the plot
    % creating a linear fit
    lm_h = fitlm(hierarchy_parc,ts_mean);

    figure
    plot(lm_h)
    xlabel('cortical hierarchy')
    ylabel('feature value')
    Lines = findobj(gca, 'Type','Line');
    CI = findobj('Type','line','DisplayName','Confidence bounds');
    CIX = CI.XData;
    CIY = CI.YData;
    lm_intercept(age) = table2array(lm_h.Coefficients(1,1));
    lm_slope(age) = table2array(lm_h.Coefficients(2,1));

    % in a nicer plot
    dist(age,:) =  lm_intercept(age)+lm_slope(age)*CIX - CIY;
    x(age,:) = [CIX CIX(end:-1:1)];
    y(age,:) = [CIY CIY(end:-1:1)+dist(age,:)*2];
end

figure('Color','w','Position',[1,1,400,500])
for age = [1:3]
    %our feature AC40: 123
    ts = mat(index_age_group(age,:),:,123);
    ts_mean = squeeze(mean(ts,1));

    % figure
    f = lm_intercept(age)+lm_slope(age)*hierarchy_parc;
    scatter(hierarchy_parc,ts_mean,90,colors{age},'filled','o','MarkerFaceAlpha',.6)
    hold on
    fill(x(age,:),y(age,:),colors{age}(1,:),'FaceAlpha',.4,'EdgeColor','none')
    h(age) = plot(hierarchy_parc,f,'LineWidth',3,'Color',colors{age}(1,:))
end

xlabel('cortical hierarchy','FontSize',16,'FontWeight','normal')
ylabel('autocorrelation, lag 133 ms','FontSize',16,'FontWeight','normal')
set(gca,'box','off')
set(gca,'FontSize',14)
hold off

exportgraphics(gcf,[filepath '\PLOTS\Fig3b_scatter_linfit_AC40_age.png'],'Resolution',500)

%% SUPPLEMENTARY: boxplots AC40 

% take a subset of the matrix mat, only at feature AC40, number: 123
mat_boxplot = squeeze(mat(:,:,123));

% make grouping arrays
group_age = [repmat(1,50,1); repmat(2.5,50,1); repmat(4,50,1)];
group_hierarchy = categorical(repelem([{'1'}, {'2'}, {'3'}], [150 150 150])');

colors = [repmat(lajolla10(4,:),50,1);repmat(lajolla10(6,:),50,1);repmat(lajolla10(8,:),50,1)];

figure('color','w', 'Position',[0,0,700,1000])
h = boxchart(repmat(group_age,3,1),[mean(mat_boxplot([1:50,151:200,301:350],idx_h_min),2);mean(mat_boxplot([1:50,151:200,301:350],idx_h_med),2);mean(mat_boxplot([1:50,151:200,301:350],idx_h_max),2)],"GroupByColor",group_hierarchy)
h(1).BoxFaceColor = lajolla10(4,:);  
h(2).BoxFaceColor = lajolla10(6,:);
h(3).BoxFaceColor = lajolla10(8,:);

hold('on');

s1 = swarmchart([repmat(0.67,50,1);repmat(1,50,1);repmat(1.33,50,1)],[mean(mat_boxplot(1:50,idx_h_min),2);mean(mat_boxplot(1:50,idx_h_med),2);mean(mat_boxplot(1:50,idx_h_max),2)],[],colors,"filled",'MarkerFaceAlpha',0.6,'XJitterWidth',0.2);
s2 = swarmchart([repmat(2.17,50,1);repmat(2.5,50,1);repmat(2.83,50,1)],[mean(mat_boxplot(151:200,idx_h_min),2);mean(mat_boxplot(151:200,idx_h_med),2);mean(mat_boxplot(151:200,idx_h_max),2)],[],colors,"filled",'MarkerFaceAlpha',0.6,'XJitterWidth',0.2);
s3= swarmchart([repmat(3.67,50,1);repmat(4,50,1);repmat(4.33,50,1)],[mean(mat_boxplot(301:350,idx_h_min),2);mean(mat_boxplot(301:350,idx_h_med),2);mean(mat_boxplot(301:350,idx_h_max),2)],[],colors,"filled",'MarkerFaceAlpha',0.6,'XJitterWidth',0.2);

xticks([1,2.5,4])
set(gca,'XTickLabel',{'young','middle','old'});
set(gca,'FontSize',14)
xlabel('age group')
ylabel('autocorrelation')
hold off

exportgraphics(gcf,[filepath '\PLOTS\SUPPL_bloxplot_AC40_age.png'],'Resolution',500)

