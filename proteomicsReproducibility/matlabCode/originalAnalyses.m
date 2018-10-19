load('Day3Exp2Data.mat') %D3_R1~rq_131_sn_scaled
load('colNames.mat')
load('Day3Exp1Data_TOUSE.mat') %D3_R2~rq_128c_sn_scaled
load('colNamesCombined.mat')

% logData = log10(Day3Exp2Data);
% logData(logData == inf) = NaN;
% logData(logData == -inf) = NaN;

% j=1
% Day3Exp2Data(Day3Exp2Data <= 1) = NaN;
% for i = 1:2:size(Day3Exp2Data,2)
%     
%     figure;
%     [r2,pvalue] = corrplot(Day3Exp2Data(:,i:i+1),'varNames',colNames(i:i+1));
% 
%     hfig = gcf;
%     haxes = findobj(hfig, 'Type', 'Axes');
%     set(haxes,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
%     set(0,'defaultfigurecolor',[1 1 1])
% 
%     filenames = {'dmso2.svg','erl2.svg','lap2.svg','soraf2.svg','sunit2.svg'};
%     print(filenames{j}, '-Painters', '-dsvg','-r600')
%     j=j+1;
% end
% 
% 
% j=1
% Day3Exp1DataTOUSE(Day3Exp1DataTOUSE <= 1) = NaN;
% for i = 1:2:size(Day3Exp1DataTOUSE,2)
%     
%     figure;
%     [r2,pvalue] = corrplot(Day3Exp1DataTOUSE(:,i:i+1),'varNames',colNames(i:i+1));
% 
%     hfig = gcf;
%     haxes = findobj(hfig, 'Type', 'Axes');
%     set(haxes,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
%     set(0,'defaultfigurecolor',[1 1 1])
% 
%     filenames = {'dmso1.svg','erl1.svg','lap1.svg','soraf1.svg','sunit1.svg'};
%     print(filenames{j}, '-Painters', '-dsvg','-r600')
%     j=j+1;
% end


% 
% for i = 1:4:size(combinedData,2)
%     
%     figure;
%     [r2,pvalue] = corrplot(combinedData(:,i:i+3),'varNames',colNamesCombined(i:i+3));
% 
%     hfig = gcf;
%     haxes = findobj(hfig, 'Type', 'Axes');
%     set(haxes,'Fontname','Arial','FontWeight','Bold','Fontsize',6)
%     set(0,'defaultfigurecolor',[1 1 1])
% 
%     filenames = {"dmso.svg","erl.svg","lap.svg","soraf.svg","sunit.svg"};
%     print(filenames{j}, '-Painters', '-dsvg','-r600')
%     j=j+1;
% end
% 
% 

%%%%%%%%%%%%%%%%%%%5
%PCA
load('Day3R1.mat') %D3_R1~rq_131_sn_scaled
load('Day3R2.mat') %D3_R2~rq_128c_sn_scaled
load('colNamesCombined_drugs.mat')


% j=1
% % combinedData = [];
% for i = 1:2:size(Day3Exp2Data,2)
%    combinedData(:,j:j+1) = Day3Exp2Data(:,i:i+1);
%    combinedData(:,j+2:j+3) = Day3Exp1DataTOUSE(:,i:i+1);
%    j =j+4;
% end

combinedData = [Day3R1 Day3R2];    
% 
% combinedData(combinedData <= 1) = 0;
% combinedData(isnan(combinedData))=0;


combinedData(combinedData <= 1) = NaN;
j=1
for row = 1:size(combinedData,1)
    if all(~isnan(combinedData(row,:)))
        combinedData2(j,:) = combinedData(row,:);
        j=j+1;
    end
end
combinedData = combinedData2;



[day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(combinedData','numcomponents',2,'Rows','complete');
% [day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(normData','numcomponents',2,'Rows','pairwise','Centered',0,'algorithm','eig');

figure('Units','centimeters', 'Position', [100, 100, 25.325, 23.031]);
% pc1 = day1coeff(:,1);
% pc2 = day1coeff(:,2);
pc1 = day1score(:,1);
pc2 = day1score(:,2);
for i = 1:2:19
    scatter(pc1(i:i+1),pc2(i:i+1),100,'filled')
    labelpoints(pc1(i:i+1),pc2(i:i+1),colNamesCombined_drugs(i:i+1))
    hold on
end

tree = linkage(combinedData','average');
figure;
dendrogram(tree,'labels',shortColNames_2)
set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',12)

tree = linkage(combinedData','average','correlation');
figure;
dendrogram(tree,'labels',shortColNames)
set(gca,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
