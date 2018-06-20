clear;clc
set(0,'DefaultFigureVisible','on');
set(0,'defaultfigurecolor',[1 1 1])


load('data/sharons_data/erl_de.mat')
load('data/sharons_data/lap_de.mat')
load('data/sharons_data/sun_de.mat')
load('data/sharons_data/sor_de.mat')
load('data/sharons_data/ensembleids_de.mat')

load('agg_time_index.mat')
load('agg_dose_index.mat')
load('index_dose_8.mat')
load('index_dose_11.mat')
load('index_time_15.mat')

full_dose_index = [agg_dose_index(1:140);index_dose_8;agg_dose_index(141:180);index_dose_11;agg_dose_index(181:end)];
full_time_index = [agg_time_index(1:260);index_time_15;agg_time_index(261:end)];

time_index = [2:5];
dose_index = [1,3,6];
time_group = [erl_de(:,time_index) ,lap_de(:,time_index),sor_de(:,time_index),sun_de(:,time_index)];
dose_group = [erl_de(:,dose_index) ,lap_de(:,dose_index),sor_de(:,dose_index),sun_de(:,dose_index)];

X_time = time_group;
X_dose = dose_group;

% time 8,5,14,12,16,10,11,2, the rest 
% dose 2,3,16,15,6,1,7,4, the rest

agg_time = X_time(full_time_index,:)';
agg_dose = X_dose(full_dose_index,:)';



% agg_time = (X_time(full_time_index,:))';
% agg_dose = (X_dose(full_dose_index,:))';
agg_dose(:,320) = 0;
% agg_time(:,292) remove?
% agg_time = [agg_time(:,1:291) agg_time(:,293:end)];

agg_dose_old = agg_dose;

newtimeorder = [8,5,14,12,16,10,11,2];
for i = 1:length(newtimeorder)
    newtimeorder_geneindex(i) = (newtimeorder(i)-1)*20+1;
end

j = 1;
k=1;

for fuck = 1:length(newtimeorder_geneindex)
    for i = 1:20:size(agg_time,2)
        if i==newtimeorder_geneindex(fuck)
            agg_time_new(:,j:j+19) = agg_time(:,i:i+19);
            j = j+20;
%         else
%             agg_time_rest(:,k:k+19) = agg_time(:,i:i+19);
%             k = k+20;
        end
    end
end

for i = 1:20:size(agg_time,2)
    if ~any(i==newtimeorder_geneindex)
        agg_time_rest(:,k:k+19) = agg_time(:,i:i+19);
        k = k+20;
    end
end
    
agg_time = 0;
agg_time = [agg_time_new agg_time_rest];



newdoseorder = [2,3,16,15,6,1,7,4];
for i = 1:length(newdoseorder)
    newdoseorder_geneindex(i) = (newdoseorder(i)-1)*20+1;
end

j = 1;
k=1;
for fuck = 1:length(newdoseorder_geneindex)
    for i = 1:20:size(agg_dose,2)
        if i==newdoseorder_geneindex(fuck)
            agg_dose_new(:,j:j+19) = agg_dose(:,i:i+19);
            j = j+20;
%         else
%             agg_dose_rest(:,k:k+19) = agg_dose(:,i:i+19);
%             k = k+20;
        end
    end
end

for i = 1:20:size(agg_dose,2)
    if ~any(i==newdoseorder_geneindex)
        agg_dose_rest(:,k:k+19) = agg_dose(:,i:i+19);
        k = k+20;
    end
end

agg_dose = 0;
agg_dose = [agg_dose_new agg_dose_rest];


% agg_dose = flipud(agg_dose);
% agg_time = flipud(agg_time);
% 
agg_dose = fliplr(agg_dose);
agg_time = fliplr(agg_time);


% agg_dose = rot180(agg_dose);
% agg_time = rot180(agg_time);






% newmap = [69,117,180;116,173,209;171,217,233;204,233,242;239,233,195;253,205,126;250,157,89;240,101,6;215,48,39];
% newmap = newmap./255;

%Slight tweak, copied from final cross plot code
newmap = [69,117,180;88,140,191;128,183,214;174,218,233;239, 233, 195;249,151,86;238,98,62;226,73,50;215,48,39];
newmap = newmap./255;
    
    
    
% figure('Position', [100, 100, 1500, 1500]);
figure('Units','centimeters', 'Position', [100, 100, 26.325, 23.031]);
agg_time(end+1,:) = 0; %Get around weird bug of cutting off last col
agg_time(:,end+1) = 0;
h=surf(agg_time');
set(h,'edgecolor','none')
colormap(newmap)
view(2)
% caxis([-3 3]);
caxis([-5.265 5.265]);
h=colorbar;
set(h,'fontsize',22,'Fontname','Arial');

    %         xlabel('Conditions, 3.16uM','Fontsize',22,'Fontname','Arial','FontWeight','Bold')
    %         xticklabels = {'Erl 6hr', 'Erl 24hr', 'Erl 72hr', 'Erl 168hr', 'Lap 6hr', 'Lap 24hr', 'Lap 72hr', 'Lap 168hr', 'Sor 6hr', 'Sor 24hr', 'Sor 72hr', 'Sor 168hr', 'Sun 6hr', 'Sun 24hr', 'Sun 72hr', 'Sun 168hr'};
    %         set(gca,'XLim',[1 17],'XTick',1.5:16.5,'XTickLabel',xticklabels,'YLim',[1 size(agg_time,2)]);
    set(gca,'XLim',[1 17],'XTickLabel','','YLim',[1 size(agg_time,2)],'YTickLabel','');
    set(gca,'Fontsize',26,'Fontname','Arial')
%     ylabel('Genes','Fontsize',34,'Fontname','Arial','FontWeight','Bold')
filename = 'aggregate_time_clusters.svg';
print(filename, '-Painters', '-dsvg','-r600')



% figure('Position', [100, 100, 1500, 1500]);
figure('Units','centimeters', 'Position', [100, 100, 26.325, 23.031]);
agg_dose(end+1,:) = 0; %Get around weird bug of cutting off last col
agg_dose(:,end+1) = 0;
h=surf(agg_dose');
set(h,'edgecolor','none')
colormap(newmap)
view(2)
% caxis([-3 3]);
caxis([-5.265 5.265]);
h=colorbar;
set(h,'fontsize',22,'Fontname','Arial');


    %         xlabel('Conditions, 24hr','Fontsize',22,'Fontname','Arial','FontWeight','Bold')
    %         xticklabels={'Erl 1uM', 'Erl 3uM', 'Erl 10 uM', 'Lap 1uM', 'Lap 3uM', 'Lap 10uM', 'Sor 1uM', 'Sor 3uM', 'Sor 10 uM', 'Sun 1uM', 'Sun 3uM', 'Sun 10uM'};
    %         set(gca,'XLim',[1 13],'XTick',1.5:12.5,'XTickLabel',xticklabels,'YLim',[1 size(agg_dose,2)]);
    set(gca,'Fontsize',26,'Fontname','Arial')
    set(gca,'XLim',[1 13],'XTickLabel','','YLim',[1 size(agg_dose,2)],'YTickLabel','');
%     ylabel('Genes','Fontsize',34,'Fontname','Arial','FontWeight','Bold')
filename = 'aggregate_dose_clusters.svg';
print(filename, '-Painters', '-dsvg','-r600')


 time_clusters = 0;
 time_labels = [16,14,13,12,11,10,9,8,7,6,5,4,3,2,1];
 for i = 1:length(time_labels)
     time_clusters = [time_clusters;ones(20,1).*time_labels(i)];
 end
 time_clusters = time_clusters(2:end);
 
 dose_clusters = 0;
 dose_labels = [16,15,14,13,12,10,9,7,6,5,4,3,2,1];
 for i = 1:length(dose_labels)
     dose_clusters = [dose_clusters;ones(20,1).*dose_labels(i)];
 end
 dose_clusters = dose_clusters(2:end);
 
 
 
timelabels = {'Cluster ID,','Erl 6hr,', 'Erl 24hr,', 'Erl 72hr,', 'Erl 168hr,', 'Lap 6hr,', 'Lap 24hr,', 'Lap 72hr,', 'Lap 168hr,', 'Sor 6hr,', 'Sor 24hr,', 'Sor 72hr,', 'Sor 168hr,', 'Sun 6hr,', 'Sun 24hr,', 'Sun 72hr,', 'Sun 168hr,'};
doselabels={'Cluster ID,','Erl 1uM,', 'Erl 3uM,', 'Erl 10 uM,', 'Lap 1uM,', 'Lap 3uM,', 'Lap 10uM,', 'Sor 1uM,', 'Sor 3uM,', 'Sor 10 uM,', 'Sun 1uM,', 'Sun 3uM,', 'Sun 10uM,'};

% filename1 = 'agg_time.csv';
% filename2 = 'agg_dose.csv';
% 
% fid1 = fopen(filename1,'w');
% fprintf(fid1, '%s', timelabels{1:end-1});
% fprintf(fid1, '%s\n',timelabels{end});
% fclose(fid1);
% 
% fid2 = fopen(filename2,'w');
% fprintf(fid2, '%s', doselabels{1:end-1});
% fprintf(fid2, '%s\n',doselabels{end});
% fclose(fid2);
% % 
% time_for_file = [time_clusters(1:end-1) flipud(agg_time(1:end-1,1:end-1))'];
% dose_for_file = [dose_clusters(1:end-1) flipud(agg_dose(1:end-1,1:end-1))'];
% 
%  dlmwrite(filename1,time_for_file,'-append');
%  dlmwrite(filename2,dose_for_file,'-append');

 
