set(0,'DefaultFigureVisible','on');
set(0,'defaultfigurecolor',[1 1 1])

%Load datasets
load('differentialExpressionData/erl_de.mat')
load('differentialExpressionData/lap_de.mat')
load('differentialExpressionData/sun_de.mat')
load('differentialExpressionData/sor_de.mat')
load('differentialExpressionData/ensembleids_de.mat')

%Load gene indices for each cluster
load('timeClustersIndex.mat')   
load('doseClustersIndex.mat')

%Add indices from clusters with no goseq results 
load('index_dose_8.mat')
load('index_dose_11.mat')
load('index_time_15.mat')
dose_index = [doseClustersIndex(1:140);index_dose_8;doseClustersIndex(141:180);index_dose_11;doseClustersIndex(181:end)];
time_index = [timeClustersIndex(1:260);index_time_15;timeClustersIndex(261:end)];

time_columns = [2:5];
dose_columns = [1,3,6];
time_group = [erl_de(:,time_columns) ,lap_de(:,time_columns),sor_de(:,time_columns),sun_de(:,time_columns)];
dose_group = [erl_de(:,dose_columns) ,lap_de(:,dose_columns),sor_de(:,dose_columns),sun_de(:,dose_columns)];

X_time = time_group;
X_dose = dose_group;

timeData = X_time(time_index,:)';
doseData = X_dose(dose_index,:)';



%reorder incorrectly ordered clusters
newtimeorder = [8,5,14,12,16,10,11,2];
for i = 1:length(newtimeorder)
    newtimeorder_geneindex(i) = (newtimeorder(i)-1)*20+1;
end

j = 1;
k=1;

for gene = 1:length(newtimeorder_geneindex)
    for i = 1:20:size(timeData,2)
        if i==newtimeorder_geneindex(gene)
            agg_time_new(:,j:j+19) = timeData(:,i:i+19);
            j = j+20;
        end
    end
end

for i = 1:20:size(timeData,2)
    if ~any(i==newtimeorder_geneindex)
        agg_time_rest(:,k:k+19) = timeData(:,i:i+19);
        k = k+20;
    end
end
    
timeData = 0;
timeData = [agg_time_new agg_time_rest];



newdoseorder = [2,3,16,15,6,1,7,4];
for i = 1:length(newdoseorder)
    newdoseorder_geneindex(i) = (newdoseorder(i)-1)*20+1;
end

j = 1;
k=1;
for gene = 1:length(newdoseorder_geneindex)
    for i = 1:20:size(doseData,2)
        if i==newdoseorder_geneindex(gene)
            agg_dose_new(:,j:j+19) = doseData(:,i:i+19);
            j = j+20;
        end
    end
end

for i = 1:20:size(doseData,2)
    if ~any(i==newdoseorder_geneindex)
        agg_dose_rest(:,k:k+19) = doseData(:,i:i+19);
        k = k+20;
    end
end

doseData = 0;
doseData = [agg_dose_new agg_dose_rest];


doseData = fliplr(doseData);
timeData = fliplr(timeData);


%Slight tweak, copied from final cross plot code
newmap = [69,117,180;88,140,191;128,183,214;174,218,233;239, 233, 195;249,151,86;238,98,62;226,73,50;215,48,39];
newmap = newmap./255;
  
%Plot time clusters    
figure('Units','centimeters', 'Position', [100, 100, 26.325, 33.031]);
timeData(end+1,:) = 0; %Get around weird bug of cutting off last row/col
timeData(:,end+1) = 0;
h=surf(timeData');
set(h,'edgecolor','none')
colormap(newmap)
view(2)
caxis([-5.265 5.265]);
h=colorbar;
set(h,'fontsize',22,'Fontname','Arial');
set(gca,'XLim',[1 17],'XTickLabel','','YLim',[1 size(timeData,2)],'YTickLabel','');
set(gca,'Fontsize',26,'Fontname','Arial')
% filename = 'time_clusters.svg';
% print(filename, '-Painters', '-dsvg','-r600')

%Plot Dose clusters
figure('Units','centimeters', 'Position', [100, 100, 26.325, 23.031]);
doseData(end+1,:) = 0; %Get around weird bug of cutting off last col
doseData(:,end+1) = 0;
h=surf(doseData');
set(h,'edgecolor','none')
colormap(newmap)
view(2)
caxis([-5.265 5.265]);
h=colorbar;
set(h,'fontsize',22,'Fontname','Arial');
set(gca,'Fontsize',26,'Fontname','Arial')
set(gca,'XLim',[1 13],'XTickLabel','','YLim',[1 size(doseData,2)],'YTickLabel','');
% filename = 'dose_clusters.svg';
% print(filename, '-Painters', '-dsvg','-r600')


