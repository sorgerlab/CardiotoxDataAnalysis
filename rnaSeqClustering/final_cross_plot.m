%Pull up eids for cluster 1 for time and dose
%look at aggregate lists
%are expression levels the same?

clear;clc
set(0,'DefaultFigureVisible','on');
set(0,'defaultfigurecolor',[1 1 1])

%load datasets
%log2 fold change
load('../data/sharons_data/erl_de.mat')
load('../data/sharons_data/lap_de.mat')
load('../data/sharons_data/sun_de.mat')
load('../data/sharons_data/sor_de.mat')
load('../data/sharons_data/ensembleids_de.mat')

%Combine datasets into single array
%Contains expression information for all differentially expressed genes,
%all condiditions
X = [erl_de,lap_de,sor_de,sun_de];

%load indices for ensemble ids for genes of interest
%These need to be generated for every dataset being plotted
load('aggregate_time_index.mat')    %These may be wrong
load('aggregate_dose_index.mat')
load('allTimeIndex.mat')    %These may be wrong
load('allDoseIndex.mat')


%Cluster indices to plot
%Already have clusters identified from previous g-means run
clusters_dose = [2,3,5,16];
clusters_time = [4,12,16];
clusters = [2,3,5,16,4,12,16];

%Pull genes from all cluster, using saved index
%Separate out dose and time conditions
agg_time2 = X(agg_time_index,:)';
agg_dose2 = X(agg_dose_index,:)';
agg_time = X(allTimeIndex,:)';
agg_dose = X(allDoseIndex,:)';
agg_dose(:,280) = 0;
agg_dose2(:,280) = 0;
% etrans = X(etransIndex,:)';

%Create cell aray with one array for each dose-specified cluster,
%containing top 20 genes at all 24 conditions
dose_cell = cell(1,16);
j = 1;
for i = 1:16 %missing 8,11
    if i == 8 || i == 11 %don't have goseq gene lists for clusters 8 and 11
        dose_cell{i} = zeros(24,20);    %Fill with zeros for now
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
    if i == 15 %don't have goseq gene lists for cluster 15
        time_cell{i} = zeros(24,20);  %Fill with zeros for now
    else
        time_cell{i} = agg_time(:,j:j+19);
        j = j+20;
    end
end


j=1;
%loop through 7 identified clusters of interest
for i = 1:length(clusters)
    %first four clusters from dose, last three from from time
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
    
    %For plotting
    %Create new custom heatmap (blue - beige - orange, 9 colors)
    %     newmap = [69,117,180;116,173,209;171,217,233;204,233,242;239,233,195;253,205,126;250,157,89;240,101,6;215,48,39];
    %     newmap = [69,117,180;171,217,233;204,233,242;239,233,195;253,205,126;250,157,89;215,48,39];
%     newmap = [69,117,180;101 156 200;138,190,218;176,219,234;239, 233, 195;253,213,134;253,178,101;247,134,78;235,90,58];
%     newmap = newmap./255;
% 
%         
%         newmap = 
%             [69,117,180
%             ;116,173,209;
%             171,217,233;
%             204,233,242;
%             239,233,195;
%             253,205,126
%             ;250,157,89;
%             240,101,6
%             ;215,48,39];

%          [69,117,180
%         101 156 200
%         138,190,218
%         176,219,234
%         213,237,245 %Probably drop this
%         239, 233, 195
%         253,213,134
%         253,178,101
%         247,134,78
%         235,90,58


%Back to R script for colors:
%R hex values:
% "#4575B4" "#588CBF" "#6BA3CB" "#80B7D6" "#97C9E0" "#AEDAE9" "#C3E5F0" "#D9EFF6" "#E8EDD9" "#F5E5AE" "#FDDA8A" "#FDC577" "#FDB063" "#F99756" "#F67C4A" "#EE623E" "#E24932" "#D73027"
% 69,117,180
% 88,140,191
% 107,163,203
% 128,183,214
% 151,201,224
% 174,218,233
% 195,229,240
% 217,239,246
% 
% 232,237,217
% 
% 245,229,174
% 253,218,138
% 253,197,119
% 253,176,99
% 249,151,86
% 246,124,74
% 238,98,62
% 226,73,50
% 215,48,39


% 69,117,180
% 88,140,191
% % % 107,163,203
% 128,183,214
% % % 151,201,224
% 174,218,233

% 232,237,217 - 239, 233, 195;
% 
% 253,197,119
% % % 253,176,99
% 249,151,86
% % % 246,124,74
% 238,98,62
% 226,73,50
% 215,48,39

    newmap = [69,117,180;88,140,191;128,183,214;174,218,233;239, 233, 195;249,151,86;238,98,62;226,73,50;215,48,39];
    newmap = newmap./255;



    %Reduce default heatmap (blue - green - yellow) to 9 colors
    %Can call custom or default colormap as desired when plotting
    colormap default
    cmap = colormap;
    defmap = cmap(1:7:64,:);
    defmap = [defmap(1:5,:);defmap(7:end,:)]; %remove 6
    
    %%%%%%%%%%%%%%%%%
    %%% Erlotinib %%%
    %%%%%%%%%%%%%%%%%
    
    %This can be left alone if number of genes (20) and conditions (6)
    %remains consistent. Much customize if these change 
    
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
   
   %Create a subplot for each cross with 4 columns and 7 rows (7 clusters x 4 drugx = 28 subplots)
   %This is figure-specific and can be altered for a different layout
   %Removing subplot calls and instead using 'figure' call for every plot
   %will generate a separate plot for each cross
   
   subplot(7,4,j);
   
   %Plot cross
   %Plotting command
   h=surf(erl_cross_array');
   %remove black mesh from plot
   set(h,'edgecolor','none')
   %remove x and y axis
   axis off
   %choose which colormap to use
%    colormap(defmap)
   colormap(newmap)
   %convert 3d plot to 2d
   view(2)
   %set color scale
   caxis([-4 4]); 
   grid off
   
   %For first cluster (first row) add name of drug as title
%    if i == 1
%        title('Erlotinib','Fontname','Arial','FontWeight','Bold','Fontsize',7)
%    end
   
   %Plot axes labels
   %Only plot y labels on first column
   %Only plot x labels on bottom row
   xticklabels={'6 h','24 h','72 h','168 h'};
   yticklabels={'10 \muM','3 \muM','1 \muM'};
   if any(j==25:28)
       set(gca,'XLim',[1 5],'XTick',1.5:2:4.5,'XTickLabel',xticklabels,'Fontname','Arial','FontWeight','Bold','Fontsize',7)
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
    axis off
    colormap(newmap) %Can choose default or custom colormap
%     colormap(defmap)
    view(2)
    caxis([-5.265 5.265]); 
    grid off
    
    %For first cluster (first row) add name of drug as title
%     if i == 1
%         title('Lapatinib','Fontname','Arial','FontWeight','Bold','Fontsize',7)
%     end
    
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
    axis off
    colormap(newmap) %Choose default or new colormap
%     colormap(defmap)
    view(2)
    caxis([-5.265 5.265]); 
    grid off
    
    %For first cluster (first row) add name of drug as title
%     if i == 1
%         title('Sorafenib','Fontname','Arial','FontWeight','Bold','Fontsize',7)
%     end

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
    axis off
    colormap(newmap) %Choose default or custom colormap
%     colormap(defmap)
    view(2)
    caxis([-5.265 5.265]); 
    grid off
    
    %For first cluster (first row) add name of drug as title
%     if i == 1
%         title('Sunitinib','Fontname','Arial','FontWeight','Bold','Fontsize',7)
%     end
    
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
    j=j+1
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
%Uncomment print command to write to file, will overwrite if there's
%another file already present with same name
%Can change file format from svg to png, etc
fn1 = 'testSVG_seven_clusters_final_default_colormap.svg';
% print(fn1, '-opengl', '-dsvg','-r600')
print(fn1, '-Painters', '-dsvg','-r600')

    