clear;clc
set(0,'DefaultFigureVisible','on');
set(0,'defaultfigurecolor',[1 1 1])

%load datasets
load('data/sharons_data/erl_de.mat')
load('data/sharons_data/lap_de.mat')
load('data/sharons_data/sun_de.mat')
load('data/sharons_data/sor_de.mat')
load('data/sharons_data/ensembleids_de.mat')

%load indices for ensemble ids for genes of interest
load('agg_time_index.mat')
load('agg_dose_index.mat')

%Cluster indices to plot
clusters_dose = [2,3,5,16];
clusters_time = [4,12,16];
clusters = [2,3,5,16,4,12,16];
%25x241, 24 conditions, 240 genes (12 clusters), row and col of zeros

%GO Term names for each cluster to plot
cluster_names = {'Sterol Biosynthetic Process/Cholesterol Biosynthetic Process','Small Molecule Biosynthetic Process/Cellular Amino Acid Metabolic Process','Response to Endoplasmic Reticulum Stress','Sarcomere/Contractile Fiber Part','Extracellular Space','Electron Transport Chain','Mitotic Cell Cycle/Mitotoic CelL Cycle Process'};

%Combine datasets into single array
X = [erl_de,lap_de,sor_de,sun_de];

%Pull genes from each cluster, using saved index
%Separate out dose and time conditions
agg_time = X(agg_time_index,:)';
agg_dose = X(agg_dose_index,:)';
agg_dose(:,280) = 0;

%Create new custom heatmap (blue - beige - orange, 9 colors)
newmap = [69,117,180;116,173,209;171,217,233;204,233,242;239,233,195;253,205,126;250,157,89;240,101,6;215,48,39];
newmap = newmap./255;

%Reduce default heatmap (blue - green - yellow) to 9 colors
%Can call custom or default colormap as desired when plotting
colormap default
cmap = colormap;
defmap = cmap(1:7:64,:);
defmap = [defmap(1:5,:);defmap(7:end,:)]; %remove 6

%Create cell aray with one array for each dose-specified cluster,
%containing top 20 genes at all 24 conditions
dose_cell = cell(1,16);
j = 1;
for i = 1:16 %missing 8,11
    if i == 8 || i == 11
        dose_cell{i} = zeros(24,20);
    else
        dose_cell{i} = agg_dose(:,j:j+19);
        j = j+20;
    end
end


%Create cell aray with one array for each time-specified cluster,
%containing top 20 genes at all 24 conditions
time_cell = cell(1,16);
j = 1;
for i = 1:16 %missing 15
    if i == 15
        time_cell{i} = zeros(24,20);
    else
        time_cell{i} = agg_time(:,j:j+19);
        j = j+20;
    end
end


j=1;
%loop through 7 identified clusters of interest
for i = 1:length(clusters)
    %first four clusters from dose, last from from time
    if i < 5
        agg_cell = dose_cell;
    else
        agg_cell = time_cell;
    end
        
    %Split out data columns relevant to each drug for individual plotting
    erl = agg_cell{clusters(i)}(1:6,:);
    lap = agg_cell{clusters(i)}(7:12,:);
    sor = agg_cell{clusters(i)}(13:18,:);
    sun = agg_cell{clusters(i)}(19:24,:);

    %%%%%%%%%%%%%%%%%
    %%% Erlotinib %%%
    %%%%%%%%%%%%%%%%%
    
    %Build array for plotting drug. Inf values will appear as white,
    %matching background. Cross of actual data will be colored
    erl_cross_array = inf(4,60);
    erl_cross_array(2,1:20) = erl(6,:);
    erl_cross_array(:,21:40) = erl(2:5,:);
    erl_cross_array(2,41:60) = erl(1,:);
    erl_cross_array(3,1:20) = 0;
    erl_cross_array(3,40:60) = 0;
     
   %Create new figure object at first cluster. Will be reused for each additional cluster subplots 
   if j == 1
        figure('Units','centimeters', 'Position', [100, 100, 15.325, 13.031]);
   end
   
   %Add row and column of zeros for account for surf bug cutting off a row
   %and col
   erl_cross_array(end+1,:) = 0;
   erl_cross_array(:,end+1) = 0;
   
   %Create a subplot for each cross (7 clusters x 4 drugx = 28 subplots)
   subplot(7,4,j);
   
   %Plot cross
   h=surf(erl_cross_array');
   set(h,'edgecolor','none')
   colormap(newmap)
   view(2)
   caxis([-3 3]);
   grid off
   
   %For first cluster (first row) add name of drug as title
   if i == 1
       title('Erlotinib','Fontname','Arial','FontWeight','Bold','Fontsize',7.0)
   end
   
   %Plot axes labels
   %Only plot y labels on first column
   %Only plot x labels on bottom row
   xticklabels={'6 h','24 h','72 h','168 h'};
   yticklabels={'10 \muM','3 \muM','1 \muM'};
   if any(j==25:28)
       set(gca,'XLim',[1 5],'XTick',1.5:4.5,'XTickLabel',xticklabels,'Fontname','Arial','FontWeight','Bold','Fontsize',7)
   else
       set(gca,'XLim',[1 5],'XTick',1.5:4.5,'XTickLabel','')
   end
   
   if any(j==1:4:28)
       set(gca,'YLim',[1 61],'YTick',10:20:50,'YTickLabel',yticklabels,'Fontname','Arial','FontWeight','Bold','Fontsize',7)
   else
       set(gca,'YLim',[1 61],'YTick',10:20:50,'YTickLabel','')
   end
   j=j+1;
    
   %Save individual crosses if desired
%    fn1 = 'heatmaps/cluster_';
%    fn2 = '/erl.png';
%    filename = [fn1 num2str(clusters(i)) fn2];
%    print(filename, '-dpng','-r300')
    

    %%%%%%%%%%%%%%%%%
    %%% Lapatinib %%%
    %%%%%%%%%%%%%%%%%
    
    %Build array for plotting drug. Inf values will appear as white,
    %matching background. Cross of actual data will be colored
    lap_cross_array = inf(4,60);
    lap_cross_array(2,1:20) = lap(6,:);
    lap_cross_array(:,21:40) = lap(2:5,:);
    lap_cross_array(2,41:60) = lap(1,:);
    lap_cross_array(3,1:20) = 0;
    lap_cross_array(3,40:60) = 0;
     
    %Add row and column of zeros for account for surf bug cutting off a row
    %and col
    lap_cross_array(end+1,:) = 0;
    lap_cross_array(:,end+1) = 0;
    
    %Create a subplot for each cross (7 clusters x 4 drugx = 28 subplots)
    subplot(7,4,j);
        
    %Plot cross
    h=surf(lap_cross_array');
    set(h,'edgecolor','none')
    %     colormap(newmap) %Can choose default or custom colormap
    colormap(newmap)
    view(2)
    caxis([-3 3]);
    grid off
    
    %For first cluster (first row) add name of drug as title
    if i == 1
        title('Lapatinib','Fontname','Arial','FontWeight','Bold','Fontsize',7)
    end
    
    %Plot axes labels
    %Only plot y labels on first column
    %Only plot x labels on bottom row
    xticklabels={'6 h','24 h','72 h','168 h'};
    yticklabels={'10 \muM','3 \muM','1 \muM'};
    if any(j==25:28)
        set(gca,'XLim',[1 5],'XTick',1.5:4.5,'XTickLabel',xticklabels,'Fontname','Arial','FontWeight','Bold','Fontsize',7)
    else
        set(gca,'XLim',[1 5],'XTick',1.5:4.5,'XTickLabel','')
    end
    
    if any(j==1:4:28)
        set(gca,'YLim',[1 61],'YTick',10:20:50,'YTickLabel',yticklabels,'Fontname','Arial','FontWeight','Bold','Fontsize',7)
    else
        set(gca,'YLim',[1 61],'YTick',10:20:50,'YTickLabel','')
    end
     j=j+1;
    
     %Save individual crosses if desired
%     fn1 = 'heatmaps/cluster_';
%     fn2 = '/lap.png';
%     filename = [fn1 num2str(clusters(i)) fn2];
% %     print(filename, '-dpng','-r300');


    %%%%%%%%%%%%%%%%%
    %%% Sorafenib %%%
    %%%%%%%%%%%%%%%%%

    %Build array for plotting drug. Inf values will appear as white,
    %matching background. Cross of actual data will be colored
    sor_cross_array = inf(4,60);
    sor_cross_array(2,1:20) = sor(6,:);
    sor_cross_array(:,21:40) = sor(2:5,:);
    sor_cross_array(2,41:60) = sor(1,:);
    sor_cross_array(3,1:20) = 0;
    sor_cross_array(3,40:60) = 0;
    
    %Add row and column of zeros for account for surf bug cutting off a row
    %and col
    sor_cross_array(end+1,:) = 0;
    sor_cross_array(:,end+1) = 0;
    
    %Create a subplot for each cross (7 clusters x 4 drugx = 28 subplots)
    subplot(7,4,j);

    %Plot Cross
    h=surf(sor_cross_array');
    set(h,'edgecolor','none')
    %     colormap(newmap) %Choose default or new colormap
    colormap(newmap)
    view(2)
    caxis([-3 3]);
    grid off
    
    %For first cluster (first row) add name of drug as title
    if i == 1
        title('Sorafenib','Fontname','Arial','FontWeight','Bold','Fontsize',7)
    end

    %Plot axes labels
    %Only plot y labels on first column
    %Only plot x labels on bottom row
    xticklabels={'6 h','24 h','72 h','168 h'};
    yticklabels={'10 \muM','3 \muM','1 \muM'};
    
    if any(j==25:28)
        set(gca,'XLim',[1 5],'XTick',1.5:4.5,'XTickLabel',xticklabels,'Fontname','Arial','FontWeight','Bold','Fontsize',7)
    else
        set(gca,'XLim',[1 5],'XTick',1.5:4.5,'XTickLabel','')
    end
    
    if any(j==1:4:28)
        set(gca,'YLim',[1 61],'YTick',10:20:50,'YTickLabel',yticklabels,'Fontname','Arial','FontWeight','Bold','Fontsize',7)
    else
        set(gca,'YLim',[1 61],'YTick',10:20:50,'YTickLabel','')
    end
    j=j+1;
    
    %Save individual crosses if desired
%     fn1 = 'heatmaps/cluster_';
%     fn2 = '/sor.png';
%     filename = [fn1 num2str(clusters(i)) fn2];
%     %     print(filename, '-dpng','-r300')
    
    
    %%%%%%%%%%%%%%%%%
    %%% Sunitinib %%%
    %%%%%%%%%%%%%%%%%

    %Build array for plotting drug. Inf values will appear as white,
    %matching background. Cross of actual data will be colored
    sun_cross_array = inf(4,60);
    sun_cross_array(2,1:20) = sun(6,:);
    sun_cross_array(:,21:40) = sun(2:5,:);
    sun_cross_array(2,41:60) = sun(1,:);
    sun_cross_array(3,1:20) = 0;
    sun_cross_array(3,40:60) = 0;
     
    %Add row and column of zeros for account for surf bug cutting off a row
    %and col
    sun_cross_array(end+1,:) = 0;
    sun_cross_array(:,end+1) = 0;
    
    %Create a subplot for each cross (7 clusters x 4 drugx = 28 subplots)
    subplot(7,4,j);
 
    %Plot cross
    h=surf(sun_cross_array');
    set(h,'edgecolor','none')
    %     colormap(newmap) %Choose default or custom colormap
    colormap(newmap)
    view(2)
    caxis([-3 3]);
    grid off
    
    %For first cluster (first row) add name of drug as title
    if i == 1
        title('Sunitinib','Fontname','Arial','FontWeight','Bold','Fontsize',7)
    end
    
    %Plot axes labels
    %Only plot y labels on first column
    %Only plot x labels on bottom row
    xticklabels={'6 h','24 h','72 h','168 h'};
    yticklabels={'10 \muM','3 \muM','1 \muM'};
    
    if any(j==25:28)
        set(gca,'XLim',[1 5],'XTick',1.5:4.5,'XTickLabel',xticklabels,'Fontname','Arial','FontWeight','Bold','Fontsize',7)
    else
        set(gca,'XLim',[1 5],'XTick',1.5:4.5,'XTickLabel','')
    end
    
    if any(j==1:4:28)
        set(gca,'YLim',[1 61],'YTick',10:20:50,'YTickLabel',yticklabels,'Fontname','Arial','FontWeight','Bold','Fontsize',7)
    else
        set(gca,'YLim',[1 61],'YTick',10:20:50,'YTickLabel','')
    end
    
    %Save individual crosses if desired
%     fn1 = 'heatmaps/cluster_';
%     fn2 = '/sun.png';
%     filename = [fn1 num2str(clusters(i)) fn2];
% %     print(filename, '-dpng','-r300')

    %Position color bar for final figure
    hp4 = get(subplot(7,4,20),'Position');
    h=colorbar('Position', [hp4(1)+hp4(3)+0.015  hp4(2)-.175  0.0275  hp4(2)+hp4(3)*2.1]);
    set(h,'Fontname','Arial','FontWeight','Bold','Fontsize',7);
    
end

%Save final figure with all 28 crosses
fn1 = 'seven_clusters_final_default_colormap.svg';
% print(fn1, '-opengl', '-dsvg','-r600')

    
%Save ensemble ids for each cluster to file
%Set file name based on conditions being clustered and cluster number
% if time == 1
%     cn1 = 'ensembleids/time/genes_cluster_';
% elseif dose == 1
%     cn1 = 'ensembleids/dose/genes_cluster_';
% elseif all == 1
%     cn1 = 'ensembleids_all/genes_cluster_';
% end
% cn2 = '.csv';
% 
% %Create cell for ensemble ids. One row for each cluster
% geneids = cell(7,1);
% %Loop through clusters
% for i = 1:7
%     %Loop through labels for cluster
%     for j = 1:length(y)
%         %Add ensemble id for gene if cluster label matches cluster index
%         if y(j) == i
%            geneids{i} = [geneids{i};ensembleids_de(j)]; 
%         end
%     end
%     %Save ids to csv
%     filename = [cn1 num2str(i) cn2];
%     dlmwrite(filename,geneids{i},'precision','%f');
% end
% 
% 
% 
% 
% 
