function [loadings,scores,loadingsTKI,scoresTKI] = plotPCA()

warning off
load('doxData.mat')
load('tkiData.mat')
load('colNamesTKIs.mat')
load('colNamesTKIsDox_subset.mat')

[normData] = normalizeData(tkiData,doxData);
normTKIOnly = normData(:,1:24);


[loadings,scores,latent,tsquared,explained] = pca(normData','numcomponents',2,'Rows','pairwise','Centered',0,'algorithm','eig');
pc1 = scores(:,1);
pc2 = scores(:,2);
pc1Var = explained(1);
pc2Var = explained(2);

% colors = {[85,98,112]./255,[78,205,196]./255,[199,244,100]./255,[255,107,107]./255,[196,77,88]./255};
% colors = {[106,74,60]./255,[233,127,2]./255,[199,244,100]./255,[255,107,107]./255,[196,77,88]./255};
% colors = {[0,160,176]./255,[106,74,60]./255,[204,51,63]./255,[235,104,65]./255,[237,201,81]./255};
colors = {[106,74,60]./255,[78,175,196]./255,[199,244,100]./255,[255,107,107]./255,[196,77,88]./255};
figure('Units','centimeters', 'Position', [100, 100, 25.325, 23.031]);
set(0,'defaultfigurecolor',[1 1 1])

drug=1;
for i = 1:6:19
    scatter(pc1(i:i+5),pc2(i:i+5),150,colors{drug},'filled')
    labelpoints(pc1(i:i+5),pc2(i:i+5),colNamesTKIsDox_subset(i:i+5),'FontSize',12)
    hold on
    drug=drug+1;
end
scatter(pc1(25:end),pc2(25:end),150,'filled')
labelpoints(pc1(25:end),pc2(25:end),colNamesTKIsDox_subset(25:end),'FontSize',12)


set(0,'defaultfigurecolor',[1 1 1])
set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',16)
xlabel(['PC1 (' sprintf('%2.1f',pc1Var) '%)'])
ylabel(['PC2 (' sprintf('%2.1f',pc2Var) '%)'])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%No geo data sets

normTKIOnly = normData(:,1:24);
[loadingsTKI,scoresTKI,latentTKI,tsquaredTKI,explainedTKI] = pca(normTKIOnly','numcomponents',2,'Rows','pairwise','Centered',0,'algorithm','eig');
pc1 = scoresTKI(:,1);
pc2 = scoresTKI(:,2);
pc1Var = explainedTKI(1);
pc2Var = explainedTKI(2);

% colors = {[85,98,112]./255,[78,205,196]./255,[199,244,100]./255,[255,107,107]./255,[196,77,88]./255};
% colors = {[106,74,60]./255,[233,127,2]./255,[199,244,100]./255,[255,107,107]./255,[196,77,88]./255};
% colors = {[0,160,176]./255,[106,74,60]./255,[204,51,63]./255,[235,104,65]./255,[237,201,81]./255};
colors = {[106,74,60]./255,[78,175,196]./255,[199,244,100]./255,[255,107,107]./255,[196,77,88]./255};
figure('Units','centimeters', 'Position', [100, 100, 25.325, 23.031]);
set(0,'defaultfigurecolor',[1 1 1])

drug=1;
for i = 1:6:19
    scatter(pc1(i:i+5),pc2(i:i+5),150,colors{drug},'filled')
    labelpoints(pc1(i:i+5),pc2(i:i+5),colNamesTKIs(i:i+5),'FontSize',12)
    hold on
    drug=drug+1;
end


set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',16)
xlabel(['PC1 (' sprintf('%2.1f',pc1Var) '%)'])
ylabel(['PC2 (' sprintf('%2.1f',pc2Var) '%)'])


