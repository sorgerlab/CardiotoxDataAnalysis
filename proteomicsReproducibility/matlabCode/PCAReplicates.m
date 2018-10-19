set(0,'DefaultFigureVisible','on');

% load('day3_allReplicates.mat')
load('synapse_day3_allReplicates.mat')
load('replicate_labels.mat')

% logData = log10(Day3Exp2Data);
% logData(logData == inf) = NaN;
% logData(logData == -inf) = NaN;


%Scale data to mean 0, std 1 
day3_allReplicates(day3_allReplicates==inf)=NaN;
day3_allReplicates(day3_allReplicates==-inf)=NaN;
day3_allReplicates(day3_allReplicates==0)=NaN;
day3_allReplicates(any(isnan(day3_allReplicates), 2), :) = [];
[normData] = normalizeData(day3_allReplicates);


%run PCA on scaled data
[loadings,scores,latent,tsquared,explained] = pca(normData','numcomponents',3,'Rows','pairwise','Centered',0,'algorithm','eig');


pc1 = scores(:,1);
pc2 = scores(:,2);
pc3 = scores(:,3);

pc1Var = explained(1);
pc2Var = explained(2);
pc3Var = explained(3);

%Build custom color scheme
colors = {[106,74,60]./255,[78,175,196]./255,[199,244,100]./255,[255,107,107]./255,[196,77,88]./255};
figure('Units','centimeters', 'Position', [100, 100, 25.325, 23.031]);
set(0,'defaultfigurecolor',[1 1 1])


%Plot first two PCs 
%For loop allows for plotting one drug at a time, 
%allowing a consistent color scheme for each drug
%Selectively label data points to improve readability 
drug=1;
for i = 1:4:size(normData,2)
    scatter(pc1(i:i+3),pc2(i:i+3),150,colors{drug},'filled')
    labelpoints(pc1(i:i+3),pc2(i:i+3),replicate_labels(i:i+3),'FontSize',12)
    hold on
    drug=drug+1;
end

%Set plot font options and label
%Add variance explained by each PC to axes labels
set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',16)
xlabel(['PC1 (' sprintf('%2.1f',pc1Var) '%)'])
ylabel(['PC2 (' sprintf('%2.1f',pc2Var) '%)'])

% print('Figures/2017/Clustering/allReplicates_PC1_PC2.png', '-Painters', '-dpng','-r600')
print('Figures/2016/PCA/allReplicates_PC1_PC2.png', '-Painters', '-dpng','-r600')



figure('Units','centimeters', 'Position', [100, 100, 25.325, 23.031]);
drug=1;
for i = 1:3
    figure
    scatter(1:size(loadings,1),loadings(:,i))
end


