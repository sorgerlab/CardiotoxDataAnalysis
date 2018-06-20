% clear
% load('colNames.mat')
% load('dox1.mat')
% load('dox2.mat')
% load('dox3.mat')
% load('dox4.mat')
% load('dox5.mat')
% load('dox6.mat')
% load('RNASeq.mat')
% 
% jointData = [RNASeq dox1 dox2 dox3 dox4 dox5 dox6];
% 
% jointData(jointData==NaN)=0;
% jointData(jointData==inf)=0;
% jointData(jointData==-inf)=0;
% 
% doxLabels = {'GSE12260.txt','GSE37260.csv','GSE42177.txt','GSE81448.txt','GSE76314','GSE40289'};
% colNames (25:30) = doxLabels;
% % 
% % [day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(jointData);
% % figure('Units','centimeters', 'Position', [100, 100, 25.325, 23.031]);
% % pc1 = day1coeff(:,1);
% % pc2 = day1coeff(:,2);
% % for i = 1:6:19
% %     scatter(pc1(i:i+5),pc2(i:i+5),100,'filled')
% %     labelpoints(pc1(i:i+5),pc2(i:i+5),colNames(i:i+5))
% %     hold on
% % end
% % scatter(pc1(25:end),pc2(25:end),100,'filled')
% % labelpoints(pc1(25:end),pc2(25:end),colNames(25:end))
% % 
% % pc1Var = day1explained(1);
% % pc2Var = day1explained(2);
% % 
% % 
% % set(0,'defaultfigurecolor',[1 1 1])
% % set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
% % xlabel(['PC1 (' num2str(pc1Var) '%)'])
% % ylabel(['PC2 (' num2str(pc2Var) '%)'])
% % 
% % fn1 = 'allData.svg';
% % print(fn1, '-Painters', '-dsvg','-r600')
% 
% % xticklabels={'6 h','24 h','72 h','168 h'};
% % yticklabels={'10 \muM','3 \muM','1 \muM'};
% % set(gca,'XLim',[1 5],'XTick',1.5:2:4.5,'XTickLabel',xticklabels,'Fontname','Arial','FontWeight','Bold','Fontsize',7)
% 
% % figure
% % bar(day1explained)
% %        
% %        
% % jointData(:,26) = []; 
% % colNames(26)=[];
% % 
% % jointData(:,end+1) = 0;
% % colNames(end+1) = 'Null';
% % 
% % [day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(jointData,'numcomponents',2);
% % figure('Units','centimeters', 'Position', [100, 100, 25.325, 23.031]);
% % pc1 = day1coeff(:,1);
% % pc2 = day1coeff(:,2);
% % % pc1 = day1score(1,:);
% % % pc2 = day1score(2,:);
% % for i = 1:6:19
% %     scatter(pc1(i:i+5),pc2(i:i+5),100,'filled')
% %     labelpoints(pc1(i:i+5),pc2(i:i+5),colNames(i:i+5))
% %     hold on
% % end
% % scatter(pc1(25:end),pc2(25:end),100,'filled')
% % labelpoints(pc1(25:end),pc2(25:end),colNames(25:end))
% % 
% % pc1Var = day1explained(1);
% % pc2Var = day1explained(2);
% % 
% % 
% % set(0,'defaultfigurecolor',[1 1 1])
% % set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
% % xlabel(['PC1 (' num2str(pc1Var) '%)'])
% % ylabel(['PC2 (' num2str(pc2Var) '%)'])
% % 
% % fn1 = 'removeOutlier.svg';
% % print(fn1, '-Painters', '-dsvg','-r600')
% % 
% % 
% % 
% % % 
% % % 
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% jointData = [RNASeq dox1 dox2 dox3 dox4 dox5 dox6];
% 
% jointData(isnan(jointData))=0;
% jointData(jointData==inf)=0;
% jointData(jointData==-inf)=0;
% 
% 
% 
% doxLabels = {'GSE12260.txt','GSE37260.csv','GSE42177.txt','GSE81448.txt','GSE76314','GSE40289'};
% colNames (25:30) = doxLabels;
%        
% % jointData(:,26) = []; 
% % colNames(26)=[];
% 
% 
% jointDataReduced = [jointData(:,2) jointData(:,8) jointData(:,14) jointData(:,20) jointData(:,25:end)];
% % colNames = {'Erl','Lap','Sor','Sun','GSE12260.txt','GSE42177.txt','GSE81448.txt','GSE76314','GSE40289'};
% colNames = {'Erl','Lap','Sor','Sun','GSE12260.txt','GSE37260.csv','GSE42177.txt','GSE81448.txt','GSE76314','GSE40289'};
% 
% 
% % jointDataReduced(isnan(jointDataReduced))=0;
% jointDataReduced(jointDataReduced==inf)=NaN;
% jointDataReduced(jointDataReduced==-inf)=NaN;
% 
% % jointDataReduced(any(isnan(jointDataReduced), 2), :) = [];
% 
% 
% % 763, 814 421
% 
% % jointDataReduced(:,5) = []; 
% % colNames(5)=[];
% % jointDataReduced(:,5) = []; 
% % colNames(5)=[];
% % jointDataReduced(:,8) = []; 
% % colNames(8)=[];
% % 
% 
% 
% jointDataReduced(:,6) = []; 
% colNames(6)=[];
% 
% % jointDataReduced(:,8) = []; 
% % colNames(8)=[];
% 
% 
% % jointDataReduced(:,9) = []; 
% % colNames(9)=[];
% % jointDataReduced(:,end+1) = 0;
% % colNames{end+1} = 'Null';
% 
% 
% % [day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(jointDataReduced','Rows','pairwise','centered','off','numcomponents',2);
% % [day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(jointDataReduced','Rows','pairwise','centered','off','numcomponents',2);
% [day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(jointDataReduced','numcomponents',2,'Rows','pairwise','algorithm','eig');
% 
% figure('Units','centimeters', 'Position', [100, 100, 25.325, 23.031]);
% % pc1 = day1coeff(:,1);
% % pc2 = day1coeff(:,2);
% pc1 = day1score(:,1);
% pc2 = day1score(:,2);
% for i = 1:4
%     scatter(pc1(i),pc2(i),100,'filled')
%     labelpoints(pc1(i),pc2(i),colNames(i))
%     hold on
%  end
% scatter(pc1(5:end),pc2(5:end),100,'filled')
% labelpoints(pc1(5:end),pc2(5:end),colNames(5:end))
% 
% pc1Var = day1explained(1);
% pc2Var = day1explained(2);

% 
% set(0,'defaultfigurecolor',[1 1 1])
% set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
% xlabel(['PC1 (' num2str(pc1Var) '%)'])
% ylabel(['PC2 (' num2str(pc2Var) '%)'])
% % 
% % figure
% %  bar(day1explained)
% 
% % figure
% % for i = 1:size(jointData,2)
% %     histogram(jointData(:,i),100)
% %     hold on
% % end
% % figure
% % histogram(jointData(:,26),100)
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCALED DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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



jointDataReduced = [jointData(:,2) jointData(:,8) jointData(:,14) jointData(:,20) jointData(:,25:end)];
% colNames = {'Erl','Lap','Sor','Sun','GSE12260.txt','GSE42177.txt','GSE81448.txt','GSE76314','GSE40289'};
colNames = {'Erl','Lap','Sor','Sun','GSE12260','GSE42177.txt','GSE81448.txt','GSE40289','GSE64476','GSE97642'};

% old
% doxLabels = {'GSE12260.txt','GSE37260.csv','GSE42177.txt','GSE81448.txt','GSE76314','GSE40289'};
%GSE37260 was outlier, 76314 poorly normalized

% jointDataReduced(isnan(jointDataReduced))=0;
jointDataReduced(jointDataReduced==inf)=NaN;
jointDataReduced(jointDataReduced==-inf)=NaN;

    
% jointDataReduced(any(isnan(jointDataReduced), 2), :) = [];


%Norm data
for row = 1:size(jointDataReduced,1)
    centeredData(row,:) = jointDataReduced(row,:)-nanmean(jointDataReduced(row,:));
end

for row = 1:size(jointDataReduced,1)
    normData(row,:) = centeredData(row,:)./nanstd(jointDataReduced(row,:));
end

normData(normData==inf)=NaN;
normData(normData==-inf)=NaN;
normData(isnan(normData))=0;

for i = 1:size(normData,1)
    if ~isnan(normData(i,end))
        normData2(i,:) = normData(i,:);
    end
end
normData = normData2;

% [day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(normData','numcomponents',2,'Rows','pairwise','Centered',0,'algorithm','svd');
[day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(normData','numcomponents',2,'Rows','pairwise','Centered',0,'algorithm','eig');
% [day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(normData','numcomponents',2,'Rows','complete','Centered',0,'algorithm','svd');
% [day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(normData','numcomponents',2,'Rows','complete','Centered',0,'algorithm','eig');
% [day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(normData','numcomponents',2,'Rows','pairwise','Centered',0,'algorithm','als');
% [day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(normData','numcomponents',2,'Rows','complete','Centered',0,'algorithm','als');


figure('Units','centimeters', 'Position', [100, 100, 25.325, 23.031]);
% pc1 = day1coeff(:,1);
% pc2 = day1coeff(:,2);
pc1 = day1score(:,1);
pc2 = day1score(:,2);
for i = 1:4
    scatter(pc1(i),pc2(i),100,'filled')
    labelpoints(pc1(i),pc2(i),colNames(i))
    hold on
 end
scatter(pc1(5:end),pc2(5:end),100,'filled')
labelpoints(pc1(5:end),pc2(5:end),colNames(5:end))

pc1Var = day1explained(1);
pc2Var = day1explained(2);


set(0,'defaultfigurecolor',[1 1 1])
set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
xlabel(['PC1 (' num2str(pc1Var) '%)'])
ylabel(['PC2 (' num2str(pc2Var) '%)'])


[day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(normData','numcomponents',2,'Rows','complete','Centered',0,'algorithm','svd');
% [day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(normData','numcomponents',2,'Rows','complete','Centered',0,'algorithm','eig');
% [day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(normData','numcomponents',2,'Rows','pairwise','Centered',0,'algorithm','als');
% [day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(normData','numcomponents',2,'Rows','complete','Centered',0,'algorithm','als');


figure('Units','centimeters', 'Position', [100, 100, 25.325, 23.031]);
% pc1 = day1coeff(:,1);
% pc2 = day1coeff(:,2);
pc1 = day1score(:,1);
pc2 = day1score(:,2);
for i = 1:4
    scatter(pc1(i),pc2(i),100,'filled')
    labelpoints(pc1(i),pc2(i),colNames(i))
    hold on
 end
scatter(pc1(5:end),pc2(5:end),100,'filled')
labelpoints(pc1(5:end),pc2(5:end),colNames(5:end))

pc1Var = day1explained(1);
pc2Var = day1explained(2);


set(0,'defaultfigurecolor',[1 1 1])
set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
xlabel(['PC1 (' num2str(pc1Var) '%)'])
ylabel(['PC2 (' num2str(pc2Var) '%)'])


[day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(normData','numcomponents',2,'Rows','complete','Centered',0,'algorithm','eig');
% [day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(normData','numcomponents',2,'Rows','pairwise','Centered',0,'algorithm','als');
% [day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(normData','numcomponents',2,'Rows','complete','Centered',0,'algorithm','als');


figure('Units','centimeters', 'Position', [100, 100, 25.325, 23.031]);
% pc1 = day1coeff(:,1);
% pc2 = day1coeff(:,2);
pc1 = day1score(:,1);
pc2 = day1score(:,2);
for i = 1:4
    scatter(pc1(i),pc2(i),100,'filled')
    labelpoints(pc1(i),pc2(i),colNames(i))
    hold on
 end
scatter(pc1(5:end),pc2(5:end),100,'filled')
labelpoints(pc1(5:end),pc2(5:end),colNames(5:end))

pc1Var = day1explained(1);
pc2Var = day1explained(2);


set(0,'defaultfigurecolor',[1 1 1])
set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
xlabel(['PC1 (' num2str(pc1Var) '%)'])
ylabel(['PC2 (' num2str(pc2Var) '%)'])


[day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(normData','numcomponents',2,'Rows','pairwise','Centered',0,'algorithm','als');
% [day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(normData','numcomponents',2,'Rows','complete','Centered',0,'algorithm','als');


figure('Units','centimeters', 'Position', [100, 100, 25.325, 23.031]);
% pc1 = day1coeff(:,1);
% pc2 = day1coeff(:,2);
pc1 = day1score(:,1);
pc2 = day1score(:,2);
for i = 1:4
    scatter(pc1(i),pc2(i),100,'filled')
    labelpoints(pc1(i),pc2(i),colNames(i))
    hold on
 end
scatter(pc1(5:end),pc2(5:end),100,'filled')
labelpoints(pc1(5:end),pc2(5:end),colNames(5:end))

pc1Var = day1explained(1);
pc2Var = day1explained(2);


set(0,'defaultfigurecolor',[1 1 1])
set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
xlabel(['PC1 (' num2str(pc1Var) '%)'])
ylabel(['PC2 (' num2str(pc2Var) '%)'])



[day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(normData','numcomponents',2,'Rows','complete','Centered',0,'algorithm','als');


figure('Units','centimeters', 'Position', [100, 100, 25.325, 23.031]);
% pc1 = day1coeff(:,1);
% pc2 = day1coeff(:,2);
pc1 = day1score(:,1);
pc2 = day1score(:,2);
for i = 1:4
    scatter(pc1(i),pc2(i),100,'filled')
    labelpoints(pc1(i),pc2(i),colNames(i))
    hold on
 end
scatter(pc1(5:end),pc2(5:end),100,'filled')
labelpoints(pc1(5:end),pc2(5:end),colNames(5:end))

pc1Var = day1explained(1);
pc2Var = day1explained(2);


set(0,'defaultfigurecolor',[1 1 1])
set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
xlabel(['PC1 (' num2str(pc1Var) '%)'])
ylabel(['PC2 (' num2str(pc2Var) '%)'])

