       clear
load('colNames.mat')
load('dox1.mat')
load('dox2.mat')
load('dox3.mat')
load('dox4.mat')
load('dox5.mat')
load('dox6.mat')
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
% jointData(:,26) = []; 
% colNames(26)=[];
% 
% jointData(:,end+1) = 0;
% colNames{end+1} = 'Null';
% 
% [day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(jointData');
% figure('Units','centimeters', 'Position', [100, 100, 25.325, 23.031]);
% % pc1 = day1coeff(:,1);
% % pc2 = day1coeff(:,2);
% pc1 = day1score(:,1);
% pc2 = day1score(:,2);
% for i = 1:6:19
%     scatter(pc1(i:i+5),pc2(i:i+5),100,'filled')
%     labelpoints(pc1(i:i+5),pc2(i:i+5),colNames(i:i+5))
%     hold on
% end
% scatter(pc1(25:end),pc2(25:end),100,'filled')
% labelpoints(pc1(25:end),pc2(25:end),colNames(25:end))
% 
% pc1Var = day1explained(1);
% pc2Var = day1explained(2);
% 
% 
% set(0,'defaultfigurecolor',[1 1 1])
% set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
% xlabel(['PC1 (' num2str(pc1Var) '%)'])
% ylabel(['PC2 (' num2str(pc2Var) '%)'])
% 
% fn1 = 'removeOutlier.svg';
% % print(fn1, '-Painters', '-dsvg','-r600')
% 
% 
% 
% 
% 







jointData = [RNASeq dox1 dox2 dox3 dox4 dox5 dox6];

jointData(jointData==NaN)=100;
jointData(jointData==inf)=100;
jointData(jointData==-inf)=100;

doxLabels = {'GSE12260.txt','GSE37260.csv','GSE42177.txt','GSE81448.txt','GSE76314','GSE40289'};
colNames (25:30) = doxLabels;
       
% jointData(:,26) = []; 
% colNames(26)=[];


jointDataReduced = [jointData(:,2) jointData(:,8) jointData(:,14) jointData(:,20) jointData(:,25:end)];
% colNames = {'Erl','Lap','Sor','Sun','GSE12260.txt','GSE42177.txt','GSE81448.txt','GSE76314','GSE40289'};
colNames = {'Erl','Lap','Sor','Sun','GSE12260.txt','GSE37260.csv','GSE42177.txt','GSE81448.txt','GSE76314','GSE40289'};

jointDataReduced(:,6) = []; 
colNames(6)=[];



[day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(jointDataReduced');
figure('Units','centimeters', 'Position', [100, 100, 25.325, 23.031]);
% pc1 = day1coeff(1,:);
% pc2 = day1coeff(2,:);
pc1 = day1score(:,1);
pc2 = day1score(:,2);
% for i = 1:6:19
    scatter(pc1,pc2,100,'filled')
    labelpoints(pc1,pc2,colNames)
    hold on
% end
scatter(pc1(25:end),pc2(25:end),100,'filled')
labelpoints(pc1(25:end),pc2(25:end),colNames(25:end))

pc1Var = day1explained(1);
pc2Var = day1explained(2);


set(0,'defaultfigurecolor',[1 1 1])
set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
xlabel(['PC1 (' num2str(pc1Var) '%)'])
ylabel(['PC2 (' num2str(pc2Var) '%)'])

