function [] = runGmeans_saveAllClusters()

set(0,'DefaultFigureVisible','off');
set(0,'defaultfigurecolor',[1 1 1])

%Load datasets
load('differentialExpressionData/erl_de.mat')
load('differentialExpressionData/lap_de.mat')
load('differentialExpressionData/sun_de.mat')
load('differentialExpressionData/sor_de.mat')
load('differentialExpressionData/ensembleids_de.mat')

%Combine all datasets
%Now clustering on entirety of dataset
data = [erl_de, lap_de, sor_de, sun_de];

%Call g-means algorithm. max clusters = 20, alpha = 0.001
[L,C] = gmeans(data, 30, 0.001);

%Find number of clusters selected (max of labels)
k=max(L);
%Reoriante data for plotting
data=data';
%Set y as cluster labels
labels=L';

%Create custom colormap (blue - beige - orange, 9 colors)
newmap = [69,117,180;88,140,191;128,183,214;174,218,233;239, 233, 195;249,151,86;238,98,62;226,73,50;215,48,39];
newmap = newmap./255;

%Check for directories to save desired heatmaps and ensebmle id spreadsheets
%Create new directories if they don't exist
if ~exist('./heatmaps','dir')
    mkdir heatmaps
end

if ~exist('./ensembleids','dir')
    mkdir ensembleids
end
    
%Create a cell array for each cluster found
%Fill with data for the appropriate cluster
for i = 1:k
    indices = find(labels==i);
    clusters{i} = data(:,indices);
    clear indices
end

%Set file names for saving heatmaps
fn1 = 'heatmaps/heatmap_';
fn2 = '.png';
 
%Loop through each cluster
for i = 1:size(clusters,2)
    i
   
    %For each cluster, create new figure object of appropriate size
    figure('Units','centimeters', 'Position', [100, 100, 15.325, 13.031]);
    %Add row and column of zeros for account for surf bug cutting off a row
    clusters{i}(end+1,:) = 0;
    clusters{i}(:,end+1) = 0;
    %Plot cluster
    h=surf(clusters{i}');
    set(h,'edgecolor','none')
    colormap(newmap)
    view(2)
    caxis([-5.265 5.265]); 
    h=colorbar;
    set(h,'fontsize',22,'Fontname','Arial'); 
    
    %Plot axes labels. Vary depending on conditions selected for clustering
        set(gca,'XLim',[1 25],'XTickLabel','','YLim',[1 size(clusters{i},2)]);
        set(gca,'Fontsize',18,'Fontname','Arial')
    
    %Save cluster heatmap to file
    filename = [fn1 num2str(i) fn2];
    print(filename, '-dpng','-r300')
    
end



%Save ensemble ids for each cluster to file
%Set file name based on conditions being clustered and cluster number
cn1 = 'ensembleids/genes_cluster_';
cn2 = '.csv';

%Create cell for ensemble ids. One row for each cluster
geneids = cell(k,1);
%Loop through clusters
for i = 1:k
    %Loop through labels for cluster
    for j = 1:length(labels)
        %Add ensemble id for gene if cluster label matches cluster index
        if labels(j) == i
           geneids{i} = [geneids{i};ensembleids_de(j)]; 
        end
    end
    %Save ids to csv
    filename = [cn1 num2str(i) cn2];
    dlmwrite(filename,geneids{i},'precision','%f');
end