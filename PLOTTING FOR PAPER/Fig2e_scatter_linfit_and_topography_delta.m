%% lin fit for theta as this has the highest lme values for the main effect of corticla hierarchy
filepath = pwd;

load([filepath '\DATA\hierarchy_parc.csv'])

% loading theta
load([filepath '\DATA\03_specData_freqbands_norm'],"delta_norm");
delta_mean = mean(zscore(delta_norm,[],2),1);

load("vik.mat")

% creating a linear fit
lm_h = fitlm(hierarchy_parc,delta_mean);

% scatter plot
figure
plot(lm_h)
xlabel('cortical hierarchy')
ylabel('feature value')
Lines = findobj(gca, 'Type','Line');
CI = findobj('Type','line','DisplayName','Confidence bounds');
CIX = CI.XData;
CIY = CI.YData;
lm_intercept = table2array(lm_h.Coefficients(1,1));
lm_slope = table2array(lm_h.Coefficients(2,1));

% in a nicer plot
dist = lm_intercept+lm_slope*CIX - CIY;
x = [CIX CIX(end:-1:1)];
y = [CIY CIY(end:-1:1)+dist*2];

% sort hierarchy_parc according to the hierarchy
[~,idx] = sort(delta_mean,'descend');
hierarchy_parc_sort = hierarchy_parc(idx);
delta_mean_sort = delta_mean(idx);

f = lm_intercept+lm_slope*hierarchy_parc;

figure('Color','w','Position',[1,1,400,600])
scatter(hierarchy_parc_sort,delta_mean_sort,90,vik(40,:),'filled','o','MarkerFaceAlpha',.5)
hold on
plot(hierarchy_parc,f,'LineWidth',3,'Color',[0.5 0.5 0.5])
fill(x,y,[0.5 0.5 0.5],'FaceAlpha',.4,'EdgeColor','none')
xlabel('cortical hierarchy','FontSize',16)
ylabel('delta power','FontSize',16)
ylim([-1.2 1.2])
set(gca,'box','off') 
set(gca,'FontSize',20)
exportgraphics(gcf,[filepath '\PLOTS\Fig2e_scatter_linfit_delta.png'],'Resolution',500)

