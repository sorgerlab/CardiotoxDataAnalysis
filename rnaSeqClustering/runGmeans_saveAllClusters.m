function [] = runGmeans_saveAllClusters()

set(0,'DefaultFigureVisible','off');
set(0,'defaultfigurecolor',[1 1 1])

%Load datasets
load('differentialExpressionData/erl_de.mat')
load('differentialExpressionData/lap_de.mat')
load('differentialExpressionData/sun_de.mat')
load('differentialExpressionData/sor_de.mat')
load('differentialExpressionData/ensembleids_de.mat')

%Select to cluster on time conditions, dose conditions, or combined
%conditions
% if strcmp(clusteringCondition,'time')
%     time = 1;
% elseif strcmp(clusteringCondition,'dose')
%     dose = 1;
% elseif strcmp(clusteringCondition,'all')
%     all = 1;
% else
%     error("Incorrect clustering condition given. Must be 'time','dose', or 'all'")
% end


%Extract time conditions (col 2:5) or dose conditions (col 1,3,6) from
%datasets
% time_index = [2:5]; %Four time points
% dose_index = [1,3,6];   %Three doses
% time_group = [erl_de(:,time_index) ,lap_de(:,time_index),sor_de(:,time_index),sun_de(:,time_index)];
% dose_group = [erl_de(:,dose_index) ,lap_de(:,dose_index),sor_de(:,dose_index),sun_de(:,dose_index)];
data = [erl_de, lap_de, sor_de, sun_de];

%Based on selected clustering variable selected, choose appropriate data
%subset
% if time == 1
%     data = time_group;
% elseif dose == 1
%     data = dose_group;
% elseif all == 1
%     data = all_group;
% end

%Call g-means algorithm. max clusters = 20, alpha = 0.001
[L,C] = gmeans(data, 30, 0.01);

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
%     mkdir heatmaps/dose
%     mkdir heatmaps/time
end

if ~exist('./ensembleids','dir')
    mkdir ensembleids
%     mkdir ensembleids/dose
%     mkdir ensembleids/time
end
    
%Create a cell array for each cluster found
%Fill with data for the appropriate cluster
for i = 1:k
    indices = find(labels==i);
    clusters{i} = data(:,indices);
    clear indices
end

%Set file names for saving heatmaps
%Vary depenind on conditions variable set earlier
% if time == 1
%     fn1 = 'heatmaps/time/heatmap_';
% elseif dose == 1
%     fn1 = 'heatmaps/dose/heatmap_';
% % elseif all == 1
% %     fn1 = 'heatmaps_all/heatmap_';
% end

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
%     if time == 1
%         set(gca,'XLim',[1 17],'XTickLabel','','YLim',[1 size(clusters{i},2)]);
%         set(gca,'Fontsize',18,'Fontname','Arial')
%         ylabel('Genes','Fontsize',34,'Fontname','Arial','FontWeight','Bold')
%     elseif dose == 1
%         set(gca,'Fontsize',18,'Fontname','Arial')
%         set(gca,'XLim',[1 13],'XTickLabel','','YLim',[1 size(clusters{i},2)]);
%         ylabel('Genes','Fontsize',34,'Fontname','Arial','FontWeight','Bold')
%     elseif all == 1
%         xlabel('Conditions','Fontsize',18,'Fontname','Arial','FontWeight','Bold')
%         xticklabels={'E1uM24hr', 'E3uM6hr', 'E3uM24hr', 'E3uM72hr', 'E3uM168hr', 'E10uM24hr','L1uM24hr', 'L3uM6hr', 'L3uM24hr', 'L3uM72hr', 'L3uM168hr', 'L10uM24hr','So1uM24hr', 'So3uM6hr', 'So3uM24hr', 'So3uM72hr', 'So3uM168hr', 'So10uM24hr','Su1uM24hr', 'Su3uM6hr', 'Su3uM24hr', 'Su3uM72hr', 'Su3uM168hr', 'Su10uM24hr'};
        set(gca,'XLim',[1 25],'XTickLabel','','YLim',[1 size(clusters{i},2)]);
        set(gca,'Fontsize',18,'Fontname','Arial')
%         ylabel('Genes','Fontsize',30,'Fontname','Arial','FontWeight','Bold')
%     end
    
    %Save cluster heatmap to file
    filename = [fn1 num2str(i) fn2];
    print(filename, '-dpng','-r300')
    
end



%Save ensemble ids for each cluster to file
%Set file name based on conditions being clustered and cluster number
% if time == 1
%     cn1 = 'ensembleids/time/genes_cluster_';
% elseif dose == 1
%     cn1 = 'ensembleids/dose/genes_cluster_';
% elseif all == 1
%     cn1 = 'ensembleids_all/genes_cluster_';
% end
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