clear;clc
set(0,'DefaultFigureVisible','off');
set(0,'defaultfigurecolor',[1 1 1])

%Load datasets
load('data/sharons_data/erl_de.mat')
load('data/sharons_data/lap_de.mat')
load('data/sharons_data/sun_de.mat')
load('data/sharons_data/sor_de.mat')
load('data/sharons_data/ensembleids_de.mat')

%Select to cluster on time conditions, dose conditions, or combined
%conditions
time = 1;
dose = 0;
all = 0;

%Extract time conditions (col 2:5) or dose conditions (col 1,3,6) from
%datasets
time_index = [2:5];
dose_index = [1,3,6];
time_group = [erl_de(:,time_index) ,lap_de(:,time_index),sor_de(:,time_index),sun_de(:,time_index)];
dose_group = [erl_de(:,dose_index) ,lap_de(:,dose_index),sor_de(:,dose_index),sun_de(:,dose_index)];
all_group = [erl_de, lap_de, sor_de, sun_de];

%Based on selected clustering variable selected, choose appropriate data
%subset
if time == 1
    X = time_group;
elseif dose == 1
    X = dose_group;
elseif all == 1
    X = all_group;
end

%Call g-means algorithm. max clusters = 20, alpha = 0.1
[L,C] = gmeans(X, 20, 0.1);
% [C] = gmeans(X, 0.9, 'pca', 'gamma');

%Find number of clusters selected (max of labels)
k=max(L);
%Reoriante data for plotting
X=X';
%Set y as cluster labels
y=L';

%Create custom colormap (blue - beige - red, 9 colors)
% newmap = [69,117,180;116,173,209;171,217,233;204,233,242;239,233,195;253,205,126;250,157,89;240,101,6;215,48,39];
newmap = [69,117,180;171,217,233;204,233,242;239,233,195;253,205,126;250,157,89;215,48,39];
newmap = newmap./255;


%Check for directories to save desired heatmaps and ensebmle id
%spreadsheets
%Create new directories if they don't exist
if ~exist('./heatmaps','dir')
    mkdir heatmaps
    mkdir heatmaps/dose
%     mkdir heatmaps_nolabels/dose_avg
    mkdir heatmaps/time
%     mkdir heatmaps_nolabels/time_avg
end

if ~exist('./ensembleids','dir')
    mkdir ensembleids
    mkdir ensembleids/dose
    mkdir ensembleids/time
end
    
%Create a cell array for each cluster found
%Fill with data for the appropriate cluster
for i = 1:k
    indices = find(y==i);
    clusters{i} = X(:,indices);
    clear indices
end

%Set file names for saving heatmaps
%Vary depenind on conditions variable set earlier
if time == 1
    fn1 = 'heatmaps/time/heatmap_';
elseif dose == 1
    fn1 = 'heatmaps/dose/heatmap_';
% elseif all == 1
%     fn1 = 'heatmaps_all/heatmap_';
end
fn2 = '.png';
 
%Loop through each cluster
for i = 1:size(clusters,2)
    i
   
    %For each cluster, create new figure object of appropriate size
%     figure('Position', [100, 100, 4500, 4500]);
    figure;
    %Add row and column of zeros for account for surf bug cutting off a row
    %and col
    clusters{i}(end+1,:) = 0;
    clusters{i}(:,end+1) = 0;
    %Plot cluster
    h=surf(clusters{i}');
    set(h,'edgecolor','none')
    colormap(newmap)
    view(2)
    caxis([-4 4]);
    h=colorbar;
    set(h,'fontsize',22,'Fontname','Arial'); 
    
    %Plot axes labels. Vary depending on conditions selected for clustering
    if time == 1
%         xlabel('Conditions, 3.16uM','Fontsize',22,'Fontname','Arial','FontWeight','Bold')
%         xticklabels = {'Erl 6hr', 'Erl 24hr', 'Erl 72hr', 'Erl 168hr', 'Lap 6hr', 'Lap 24hr', 'Lap 72hr', 'Lap 168hr', 'Sor 6hr', 'Sor 24hr', 'Sor 72hr', 'Sor 168hr', 'Sun 6hr', 'Sun 24hr', 'Sun 72hr', 'Sun 168hr'};
%         set(gca,'XLim',[1 17],'XTick',1.5:16.5,'XTickLabel',xticklabels,'YLim',[1 size(clusters{i},2)]);
        set(gca,'XLim',[1 17],'XTickLabel','','YLim',[1 size(clusters{i},2)]);
        set(gca,'Fontsize',18,'Fontname','Arial')
        ylabel('Genes','Fontsize',34,'Fontname','Arial','FontWeight','Bold')
    elseif dose == 1
%         xlabel('Conditions, 24hr','Fontsize',22,'Fontname','Arial','FontWeight','Bold')
%         xticklabels={'Erl 1uM', 'Erl 3uM', 'Erl 10 uM', 'Lap 1uM', 'Lap 3uM', 'Lap 10uM', 'Sor 1uM', 'Sor 3uM', 'Sor 10 uM', 'Sun 1uM', 'Sun 3uM', 'Sun 10uM'};
%         set(gca,'XLim',[1 13],'XTick',1.5:12.5,'XTickLabel',xticklabels,'YLim',[1 size(clusters{i},2)]);
        set(gca,'Fontsize',18,'Fontname','Arial')
        set(gca,'XLim',[1 13],'XTickLabel','','YLim',[1 size(clusters{i},2)]);
        ylabel('Genes','Fontsize',34,'Fontname','Arial','FontWeight','Bold')
    elseif all == 1
        xlabel('Conditions','Fontsize',18,'Fontname','Arial','FontWeight','Bold')
        xticklabels={'E1uM24hr', 'E3uM6hr', 'E3uM24hr', 'E3uM72hr', 'E3uM168hr', 'E10uM24hr','L1uM24hr', 'L3uM6hr', 'L3uM24hr', 'L3uM72hr', 'L3uM168hr', 'L10uM24hr','So1uM24hr', 'So3uM6hr', 'So3uM24hr', 'So3uM72hr', 'So3uM168hr', 'So10uM24hr','Su1uM24hr', 'Su3uM6hr', 'Su3uM24hr', 'Su3uM72hr', 'Su3uM168hr', 'Su10uM24hr'};
        set(gca,'XLim',[1 25],'XTick',1.5:24.5,'XTickLabel',xticklabels,'YLim',[1 size(clusters{i},2)]);
        set(gca,'Fontsize',18,'Fontname','Arial')
%         set(gca,'XLim',[1 13],'XTickLabel',xticklabels,'YLim',[1 size(clusters{i},2)]);
        ylabel('Genes','Fontsize',30,'Fontname','Arial','FontWeight','Bold')
    end
    
    %Save cluster heatmap to file
    filename = [fn1 num2str(i) fn2];
    print(filename, '-dpng','-r300')

% %Calculate and plot mean for cluster
%     clusmean = mean(clusters{i},2);
%     clusmean(:,end+1) = 0;
%     clusmean_agg(:,i) = clusmean(:,1);
% 
%     figure('Position', [100, 100, 1500, 150]);
%     h=surf(clusmean');
%     set(h,'edgecolor','none')
%     view(2)
%     caxis([-5 5]);
%     colormap(newmap)
%     colorbar;  
%     if time == 1
% %         xlabel('Conditions, 3.16uM','Fontsize',22,'Fontname','Arial','FontWeight','Bold')
% %         xticklabels = {'Erl 6hr', 'Erl 24hr', 'Erl 72hr', 'Erl 168hr', 'Lap 6hr', 'Lap 24hr', 'Lap 72hr', 'Lap 168hr', 'Sor 6hr', 'Sor 24hr', 'Sor 72hr', 'Sor 168hr', 'Sun 6hr', 'Sun 24hr', 'Sun 72hr', 'Sun 168hr'};
% %         set(gca,'XLim',[1 17],'XTick',1.5:16.5,'XTickLabel',xticklabels,'YLim',[1 size(clusters{i},2)]);
%         set(gca,'Fontsize',14,'Fontname','Arial')
%         set(gca,'XLim',[1 17],'XTickLabel','','YLim',[1 size(clusters{i},2)]);
%         ylabel('Genes','Fontsize',22,'Fontname','Arial','FontWeight','Bold')
%     else
% %         xlabel('Conditions, 24 hr','Fontsize',22,'Fontname','Arial','FontWeight','Bold')
% %         xticklabels={'Erl 1uM', 'Erl 3uM', 'Erl 10 uM', 'Lap 1uM', 'Lap 3uM', 'Lap 10uM', 'Sor 1uM', 'Sor 3uM', 'Sor 10 uM', 'Sun 1uM', 'Sun 3uM', 'Sun 10uM'};
% %         set(gca,'XLim',[1 13],'XTick',1.5:12.5,'XTickLabel',xticklabels);
%         set(gca,'XLim',[1 13],'XTickLabel','');
%         set(gca,'Fontsize',14,'Fontname','Arial')
%         ylabel('Genes','Fontsize',22,'Fontname','Arial','FontWeight','Bold')
%         yticklabels={'1'};
%         set(gca,'Ylim',[1 2],'YTick',1,'YTickLabel',yticklabels);
%     end
%     
%     %Save heatmap for each mean cluster 
%     filename = [fn3 num2str(i) fn2];
% %     print(filename, '-dpng','-r600')

    
end


% %Plot mean of all clusters in one map
% figure('Position', [100, 100, 1500, 1500]);
% clusmean_agg(:,end+1) = 0;
%     h=surf(clusmean_agg');
%     colormap(newmap)
%     set(h,'edgecolor','none')
%     view(2),[1 17],
%     caxis([-5 5]);
%     colorbar;  
%     if time == 1
% %         xlabel('Conditions, 3.16uM','Fontsize',22,'Fontname','Arial','FontWeight','Bold')
% %         xticklabels = {'Erl 6hr', 'Erl 24hr', 'Erl 72hr', 'Erl 168hr', 'Lap 6hr', 'Lap 24hr', 'Lap 72hr', 'Lap 168hr', 'Sor 6hr', 'Sor 24hr', 'Sor 72hr', 'Sor 168hr', 'Sun 6hr', 'Sun 24hr', 'Sun 72hr', 'Sun 168hr'};
% %         set(gca,'XLim',[1 17],'XTick',1.5:16.5,'XTickLabel',xticklabels);
%         set(gca,'XLim',[1 17],'XTickLabel','');
%         set(gca,'Fontsize',14,'Fontname','Arial')
% %         ylabel('Clusters','Fontsize',22,'Fontname','Arial','FontWeight','Bold')
% %         set(gca,'Ylim',[1 k+1],'YTick',1.5:k+0.5,'YTickLabel',1:k);
%         set(gca,'Ylim',[1 k+1]);
%     else
% %         xlabel('Conditions, 24 hr','Fontsize',22,'Fontname','Arial','FontWeight','Bold')
% %         xticklabels={'Erl 1uM', 'Erl 3uM', 'Erl 10 uM', 'Lap 1uM', 'Lap 3uM', 'Lap 10uM', 'Sor 1uM', 'Sor 3uM', 'Sor 10 uM', 'Sun 1uM', 'Sun 3uM', 'Sun 10uM'};
% %         set(gca,'XLim',[1 13],'XTick',1.5:12.5,'XTickLabel',xticklabels);
%         set(gca,'Fontsize',14,'Fontname','Arial')
%         set(gca,'XLim',[1 13],'XTickLabel','',[1 13]);
% %         ylabel('Clusters','Fontsize',22,'Fontname','Arial','FontWeight','Bold')
%         set(gca,'Ylim',[1 k+1],'YTick',1.5:k+0.5,'YTickLabel',1:k);
%         set(gca,'Ylim',[1 k+1]);
%     end
%     
%      
%     filename = [fn4 fn2];
%     print(filename, '-dpng','-r600')



%Save ensemble ids for each cluster to file
%Set file name based on conditions being clustered and cluster number
if time == 1
    cn1 = 'ensembleids/time/genes_cluster_';
elseif dose == 1
    cn1 = 'ensembleids/dose/genes_cluster_';
elseif all == 1
    cn1 = 'ensembleids_all/genes_cluster_';
end
cn2 = '.csv';

%Create cell for ensemble ids. One row for each cluster
geneids = cell(k,1);
%Loop through clusters
for i = 1:k
    %Loop through labels for cluster
    for j = 1:length(y)
        %Add ensemble id for gene if cluster label matches cluster index
        if y(j) == i
           geneids{i} = [geneids{i};ensembleids_de(j)]; 
        end
    end
    %Save ids to csv
    filename = [cn1 num2str(i) cn2];
    dlmwrite(filename,geneids{i},'precision','%f');
end