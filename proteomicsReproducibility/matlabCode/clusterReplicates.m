set(0,'DefaultFigureVisible','on');

% load('day3_allReplicates.mat')
load('synapse_day3_allReplicates.mat')

load('replicate_labels.mat')

day3_allReplicates(day3_allReplicates <=1) = NaN;
day3_allReplicates(any(isnan(day3_allReplicates),2),:) = [];

tree = linkage(day3_allReplicates','average');
figure;
dendrogram(tree,'Orientation','right','labels',replicate_labels)
set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
print('Figures/2017/Clustering/Single_Linkage.png', '-Painters', '-dpng','-r600')
print('Figures/2016/Clustering/Single_Linkage.png', '-Painters', '-dpng','-r600')


tree = linkage(day3_allReplicates');
figure;
dendrogram(tree,'Orientation','right','labels',replicate_labels)
set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
% print('Figures/2017/Clustering/Average_Linkage.png', '-Painters', '-dpng','-r600')
print('Figures/2016/Clustering/Average_Linkage.png', '-Painters', '-dpng','-r600')


tree = linkage(day3_allReplicates','complete');
figure;
dendrogram(tree,'Orientation','right','labels',replicate_labels)
set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
% print('Figures/2017/Clustering/Complete_Linkage.png', '-Painters', '-dpng','-r600')
print('Figures/2016/Clustering/Complete_Linkage.png', '-Painters', '-dpng','-r600')



% 
% randOrder = randperm(size(day3_allReplicates,2));
% randData = day3_allReplicates(:,randOrder);
% randLabels = replicate_labels(randOrder);
% 
% tree = linkage(randData','average');
% figure;
% dendrogram(tree,'Orientation','right','labels',randLabels)
% set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
% print('Figures/Clustering/NegativeControl_Average_Linkage.png', '-Painters', '-dpng','-r600')
% 
% tree = linkage(randData','single');
% figure;
% dendrogram(tree,'Orientation','right','labels',randLabels)
% set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
% 
% tree = linkage(randData','complete');
% figure;
% dendrogram(tree,'Orientation','right','labels',randLabels)
% set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
% 
