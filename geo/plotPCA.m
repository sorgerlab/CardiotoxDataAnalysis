

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCALED DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('doxData.mat')
load('tkiData.mat')

[normData] = normalizeData(tkiData,doxData);

% load('colNames.mat')
% colNames = colNames(1:end-1);
load('colNamesReduced.mat')
colNames = colNamesReduced;
colNames(25:end) = {'GSE12260','GSE42177','GSE81448','GSE64476','GSE40289','GSE97642'}; %overwrite generic Dox1, dox2 labels





[loadings,scores,latent,tsquared,explained] = pca(normData','numcomponents',2,'Rows','pairwise','Centered',0,'algorithm','eig');

pc1 = day1score(:,1);
pc2 = day1score(:,2);
% pc1 = day1coeff(:,1);
% pc2 = day1coeff(:,2);

% colors = {[85,98,112]./255,[78,205,196]./255,[199,244,100]./255,[255,107,107]./255,[196,77,88]./255};
% colors = {[106,74,60]./255,[233,127,2]./255,[199,244,100]./255,[255,107,107]./255,[196,77,88]./255};
% colors = {[0,160,176]./255,[106,74,60]./255,[204,51,63]./255,[235,104,65]./255,[237,201,81]./255};
colors = {[106,74,60]./255,[78,175,196]./255,[199,244,100]./255,[255,107,107]./255,[196,77,88]./255};
figure('Units','centimeters', 'Position', [100, 100, 25.325, 23.031]);

j=1;
for i = 1:6:19
    scatter(pc1(i:i+5),pc2(i:i+5),150,colors{j},'filled')
    labelpoints(pc1(i:i+5),pc2(i:i+5),colNames(i:i+5),'FontSize',12)
    hold on
    j=j+1;
end
scatter(pc1(25:end),pc2(25:end),150,'filled')
labelpoints(pc1(25:end),pc2(25:end),colNames(25:end),'FontSize',12)


pc1Var = day1explained(1);
pc2Var = day1explained(2);


set(0,'defaultfigurecolor',[1 1 1])
set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',16)
xlabel(['PC1 (' sprintf('%2.1f',pc1Var) '%)'])
ylabel(['PC2 (' sprintf('%2.1f',pc2Var) '%)'])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%No geo data sets


clear
load('allData.mat')

jointDataReduced = RNASeq; 
load('colNames.mat')
colNames = colNames(1:end-1);
% load('colNamesReduced.mat')
% colNames = colNamesReduced;
% colNames(25:end) = {'GSE12260','GSE42177','GSE81448','GSE64476','GSE40289','GSE97642'}; %overwrite generic Dox1, dox2 labels

% 
% colNames = colNames(1:end-1);
% jointDataReduced(:,end) = [];


% colNames = colNames(1:24);
% jointDataReduced(:,25:end) = [];



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


[day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(normData','numcomponents',2,'Rows','pairwise','Centered',0,'algorithm','eig');


figure('Units','centimeters', 'Position', [100, 100, 25.325, 23.031]);
% pc1 = day1coeff(:,1);
% pc2 = day1coeff(:,2);
pc1 = day1score(:,1);
pc2 = day1score(:,2);
% colors = {'g','b','o','c'};
% 
% newmap = [69,117,180;171,217,233;204,233,242;239,233,195;253,205,126;250,157,89;215,48,39];
% newmap = newmap./255;


% colors = {[85,98,112]./255,[78,205,196]./255,[199,244,100]./255,[255,107,107]./255,[196,77,88]./255};
colors = {[106,74,60]./255,[78,175,196]./255,[199,244,100]./255,[255,107,107]./255};    %,[196,77,88]./255};
% colors = {[106,74,60]./255,[233,127,2]./255,[199,244,100]./255,[255,107,107]./255,[196,77,88]./255};

% colors = {[0,160,176]./255,[106,74,60]./255,[204,51,63]./255,[235,104,65]./255,[237,201,81]./255};

j=1;
for i = 1:6:19
    scatter(pc1(i:i+5),pc2(i:i+5),150,colors{j},'filled')
%     scatter(pc1(i:i+5),pc2(i:i+5),150,'filled')
    labelpoints(pc1(i:i+5),pc2(i:i+5),colNames(i:i+5),'FontSize',12)
    hold on
    j=j+1;
end
% scatter(pc1(25:end),pc2(25:end),150,'filled')
% labelpoints(pc1(25:end),pc2(25:end),colNames(25:end),'FontSize',12)


% for i = 1:4
%     scatter(pc1(i),pc2(i),100,'filled')
%     labelpoints(pc1(i),pc2(i),colNames(i))
%     hold on
%  end
% scatter(pc1(5:end),pc2(5:end),100,'filled')
% labelpoints(pc1(5:end),pc2(5:end),colNames(5:end))

pc1Var = day1explained(1);
pc2Var = day1explained(2);


set(0,'defaultfigurecolor',[1 1 1])
set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',16)
xlabel(['PC1 (' sprintf('%2.1f',pc1Var) '%)'])
ylabel(['PC2 (' sprintf('%2.1f',pc2Var) '%)'])



% 