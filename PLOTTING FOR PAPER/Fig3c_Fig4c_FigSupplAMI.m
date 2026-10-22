filepath = pwd;

%% PLOTTING AMI OVER DIFFERENT AGES
load('vik.mat')
load('lipari10.mat')
load('lajolla10.mat')
load("batlow10.mat")

% load ami from lag 1 to 40
load([filepath '\DATA\ami_1_40.mat'])

% load t_t-stat from ami from lag 1 to 40
load([filepath '\DATA\ami_1_40_tstat.mat'])

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

%% Figure 4b 
colors = [lajolla10(4,:),.8;lajolla10(6,:),.8;lajolla10(8,:),.8];
colorbw = [.9 .9 .9; .6 .6 .6; .3 .3 .3];

for i = 1:5; ticks{i} = num2str(round(i*7.5/300*1000));end

figure('Position', [10 10 350 1000])
for age = 1:3
    ts = ami_out(index_age_group(age,:),:,:);
    ts1 = squeeze(mean(ts,1));

    % plot the autocorrelation
    subplot(4,1,age)
   
    hold('on');
    for parcel = 1:3
        ts2 = mean(ts1(eval(hierarchy_labels{parcel}),:),1);
        plot(ts2,'-','Color',colors(parcel,:),'LineWidth',3);
    end

    set(gca, 'YScale', 'log')
    box off
    xticks([0 7.5 15 22.5 30 37.5])
    xticklabels([0 ticks])
    ylim([0 1.3])
    yticks([0.01 0.1 1])
    yticklabels({'0.01', '0.1', '1'})
    set(gca,'FontSize',14)
    set(gca, 'Color', 'None');

end

subplot(4,1,4)
for cond = 1:3 

    hold('on');
    plot(allst(:,1+cond),'-','Color',colorbw(cond,:),'LineWidth',2);
    hold('off')

    ylim([-60 50])
    xlabel('lag [ms]')
    xticks([0 7.5 15 22.5 30 37.5])
    xticklabels([0 ticks])
    set(gca,'FontSize',14)
    box off
    set(gca, 'Color', 'None');
end  

exportgraphics(gcf,[filepath '\PLOTS\Fig4b_interarction_AMI_age.png'],'Resolution',500)

%% Figure 3c AMI lag 40 across age topography

colors = {repmat(batlow10(5,:),200,1),repmat(batlow10(7,:),200,1),repmat(batlow10(9,:),200,1)};

for age = [1:3]
    ts = ami_out(index_age_group(age,:),:,40);
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
    ts = ami_out(index_age_group(age,:),:,40);
    ts_mean = squeeze(mean(ts,1));

    % figure
    f = lm_intercept(age)+lm_slope(age)*hierarchy_parc;
    scatter(hierarchy_parc,ts_mean,90,colors{age},'filled','o','MarkerFaceAlpha',.6)
    hold on
    fill(x(age,:),y(age,:),colors{age}(1,:),'FaceAlpha',.4,'EdgeColor','none')
    h(age) = plot(hierarchy_parc,f,'LineWidth',3,'Color',colors{age}(1,:))
end

xlabel('cortical hierarchy','FontSize',16,'FontWeight','normal')
ylabel('auto mutual information, lag 133 ms','FontSize',16,'FontWeight','normal')
set(gca,'box','off')
set(gca,'FontSize',14)
hold off

exportgraphics(gcf,[filepath '\PLOTS\Fig3c_scatter_linfit_AMI40_age.png'],'Resolution',500)

%% Figure 3c AMI lag 40 across age topography
data2plot = squeeze(ami_out(:,:,40));

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
imwrite(Imerged,[filepath '\PLOTS\Fig3c_AMI_topo_across_age.png']);


%% SUPPLEMENTARY: boxplots AMI40 

% only at feature AMI40, number: 40
mat_boxplot = squeeze(ami_out(:,:,40));
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
ylabel('auto mutual information')
hold off

exportgraphics(gcf,[filepath '\PLOTS\SUPPL_bloxplot_AMI40_age.png'],'Resolution',500)


