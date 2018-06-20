clear;clc
set(0,'DefaultFigureVisible','on');
set(0,'defaultfigurecolor',[1 1 1])

%Load data sets
load('data/sharons_data/erl_de.mat')
load('data/sharons_data/lap_de.mat')
load('data/sharons_data/sun_de.mat')
load('data/sharons_data/sor_de.mat')
load('data/sharons_data/ensembleids_de.mat')

% load('aggregate_combined_index.mat')
time = 1;
dose = 0;

% if time == 1
    load('agg_time_index.mat')
% else
    load('agg_dose_index.mat')
% end



time_index = [2:5];
dose_index = [1,3,6];
time_group = [erl_de(:,time_index) ,lap_de(:,time_index),sor_de(:,time_index),sun_de(:,time_index)];
dose_group = [erl_de(:,dose_index) ,lap_de(:,dose_index),sor_de(:,dose_index),sun_de(:,dose_index)];

X_time = time_group;
X_dose = dose_group;
% X_all = [erl_de,lap_de,sor_de,sun_de];

agg_time = X_time(agg_time_index,:)';
agg_dose = X_dose(agg_dose_index,:)';
% agg_time(:,292) remove?
agg_time = [agg_time(:,1:291) agg_time(:,293:end)];
% agg_all = X_all(agg_all_index,:)';

newmap = [69,117,180;116,173,209;171,217,233;204,233,242;239,233,195;253,205,126;250,157,89;240,101,6;215,48,39];
newmap = newmap./255;

% figure('Position', [100, 100, 1500, 1500]);
figure;
% agg_all(end+1,:) = 0; %Get around weird bug of cutting off last col
% agg_all(:,end+1) = 0;
% agg_all(5,1:50) = inf;
if time == 1
    agg = agg_time;
    agg(end+1,:) = 0;
    agg(:,end+1) = 0;
end


h=surf(agg');
set(h,'edgecolor','none')
colormap(newmap)
view(2)
caxis([-3 3]);
h=colorbar;
set(h,'fontsize',22,'Fontname','Arial');


    %         xlabel('Conditions, 24hr','Fontsize',22,'Fontname','Arial','FontWeight','Bold')
    %         xticklabels={'Erl 1uM', 'Erl 3uM', 'Erl 10 uM', 'Lap 1uM', 'Lap 3uM', 'Lap 10uM', 'Sor 1uM', 'Sor 3uM', 'Sor 10 uM', 'Sun 1uM', 'Sun 3uM', 'Sun 10uM'};
    %         set(gca,'XLim',[1 13],'XTick',1.5:12.5,'XTickLabel',xticklabels,'YLim',[1 size(agg_dose,2)]);
    set(gca,'Fontsize',26,'Fontname','Arial')
    set(gca,'XLim',[1 26],'XTickLabel','','YLim',[1 size(agg_all,2)]);
    ylabel('Genes','Fontsize',34,'Fontname','Arial','FontWeight','Bold')
filename = 'aggregate_dose.png';
% print(filename, '-dpng','-r300')
% 
% 
%  all_clusters = 0;
%  all_labels = [15,14,13,12,11,10,9,8,6,5,3,2,1];
%  for i = 1:length(all_labels)
%      all_clusters = [all_clusters;ones(20,1).*all_labels(i)];
%  end
%  all_clusters = all_clusters(2:end);
%  
% 
%  
% timelabels = {'Cluster ID,','Erl 6hr,', 'Erl 24hr,', 'Erl 72hr,', 'Erl 168hr,', 'Lap 6hr,', 'Lap 24hr,', 'Lap 72hr,', 'Lap 168hr,', 'Sor 6hr,', 'Sor 24hr,', 'Sor 72hr,', 'Sor 168hr,', 'Sun 6hr,', 'Sun 24hr,', 'Sun 72hr,', 'Sun 168hr,'};
% doselabels={'Cluster ID,','Erl 1uM,', 'Erl 3uM,', 'Erl 10 uM,', 'Lap 1uM,', 'Lap 3uM,', 'Lap 10uM,', 'Sor 1uM,', 'Sor 3uM,', 'Sor 10 uM,', 'Sun 1uM,', 'Sun 3uM,', 'Sun 10uM,'};
% 
% filename1 = 'agg_all.csv';
% 
% 
% fid1 = fopen(filename1,'w');
% fprintf(fid1, '%s', alllabels{1:end-1});
% fprintf(fid1, '%s\n',alllabels{end});
% fclose(fid1);
% 
% 
% all_for_file = [all_clusters(1:end-1) flipud(agg_all(1:end-1,1:end-1))'];
% 
%  dlmwrite(filename1,all_for_file,'-append');
% 
%  
