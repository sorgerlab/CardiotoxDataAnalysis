%ensembleids_de - cell array in data
%Import cross_1_dose_cluster_2.csv files gives mat in correct format -
%should be easy to combine these too
%Just need to import all cross csv plots. 

%ensembleIds_20perCluster
%Split to dose/time
%Can allow for  rebuilding whole agg_dose and time index matrices, so don't
%have to change any downstream code

clear;clc;
load('ensembleids_de.mat')
load('allDoseGenes.mat')
load('allTimeGenes.mat')

allDoseIndex = zeros(length(allDoseGenes),1);
for i = 1:length(allDoseGenes)
    for j = 1:length(ensembleids_de)
        if allDoseGenes(i) == ensembleids_de{j}
            allDoseIndex(i) = j;
        end
    end
end


allTimeIndex = zeros(length(allTimeGenes),1);
for i = 1:length(allTimeGenes)
    for j = 1:length(ensembleids_de)
        if allTimeGenes(i) == ensembleids_de{j}
            allTimeIndex(i) = j;
        end
    end
end



etransIndex = zeros(length(etransportIDs),1);
for i = 1:length(etransportIDs)
    for j = 1:length(ensembleids_de)
        if etransportIDs(i) == ensembleids_de{j}
            etransIndex(i) = j;
        end
    end
end