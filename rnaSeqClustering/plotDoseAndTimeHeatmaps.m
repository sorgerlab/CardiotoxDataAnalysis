set(0,'DefaultFigureVisible','on');
set(0,'defaultfigurecolor',[1 1 1])

load('doseData.mat')
load('timeData.mat')



%Slight tweak, copied from final cross plot code
newmap = [69,117,180;88,140,191;128,183,214;174,218,233;239, 233, 195;249,151,86;238,98,62;226,73,50;215,48,39];
newmap = newmap./255;
  
%Build clustergrams
% h=clustergram(timeData','Cluster','Column','OptimalLeafOrder',true,'Colormap',newmap,'DisplayRange',5.265);
% h=clustergram(doseData','Cluster','Column','OptimalLeafOrder',true,'Colormap',newmap,'DisplayRange',5.265);
% h=clustergram(clusterMeans','Cluster','Column','OptimalLeafOrder',true,'Colormap',newmap,'DisplayRange',5.265);


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

