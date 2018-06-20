clear
load('colNames.mat')
load('dox1.mat')
load('dox2.mat')
load('dox3.mat')
load('dox4.mat')
load('dox5.mat')
load('dox6.mat')
load('RNASeq.mat')

jointData = [RNASeq dox1 dox2 dox3 dox4 dox5 dox6];

jointData(jointData==NaN)=0;
jointData(jointData==inf)=0;
jointData(jointData==-inf)=0;

doxLabels = {}
colNames (25:30) = doxLabels;

[day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(jointData);
figure;
pc1 = day1coeff(:,1);
pc2 = day1coeff(:,2);
for i = 1:6:19
    scatter(pc1(i:i+5),pc2(i:i+5),100,'filled')
    labelpoints(pc1(i:i+5),pc2(i:i+5),colNames(i:i+5))
    hold on
end
scatter(pc1(25:end),pc2(25:end),100,'filled')
labelpoints(pc1(25:end),pc2(25:end),colNames(25:end))

pc1Var = day1explained(1);
pc2Var = day1explained(2);


set(0,'defaultfigurecolor',[1 1 1])
set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
xlabel('PC1')
ylabel('PC2')

% xticklabels={'6 h','24 h','72 h','168 h'};
% yticklabels={'10 \muM','3 \muM','1 \muM'};
% set(gca,'XLim',[1 5],'XTick',1.5:2:4.5,'XTickLabel',xticklabels,'Fontname','Arial','FontWeight','Bold','Fontsize',7)
       
       
       
% dox2(:) = 0;
