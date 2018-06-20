load('colNames_Day1.mat')
load('colNames_Day1_repRemoved.mat')
load('colNames_Day3.mat')
load('day1_IBAQ_Scaled.mat')
load('day3_IBAQ_Scaled.mat')
load('day1_IBAQ_Scaled_RepRemoved.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%IBAQ Normalized Scalled to 100% - Day 1


% hfig = gcf;
% haxes = findobj(hfig, 'Type', 'Axes');
% arrayfun(@(ax) xlim(ax, [0 25]), haxes);
% arrayfun(@(ay) ylim(ay, [0 25]), haxes)


[day1coeff,day1score,day1latent,day1tsquared,day1explained] = pca(day1_IBAQ_Scaled);
figure;
pc1 = day1coeff(:,1);
pc2 = day1coeff(:,2);
for i = 1:2:9
    scatter(pc1(i:i+1),pc2(i:i+1),100,'filled')
    labelpoints(pc1(i:i+1),pc2(i:i+1),colNames_Day1(i:i+1))
    hold on
end

tree = linkage(day1_IBAQ_Scaled');
figure;
dendrogram(tree,'labels',colNames_Day1)

day1_IBAQ_Scaled(day1_IBAQ_Scaled>=30)=NaN;
[r,pvalue] = corrplot(day1_IBAQ_Scaled,'varNames',colNames_Day1,'rows','complete');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%IBAQ Normalized Scalled to 100% - Day 1 - Bad replicate removed


[day1coeff_rr,day1score_rr,day1latent_rr,day1tsquared_rr,day1explained_rr] = pca(day1_IBAQ_Scaled_RepRemoved);
figure;
pc1 = day1coeff_rr(:,1);
pc2 = day1coeff_rr(:,2);
scatter(pc1(1),pc2(1),100,'filled')
labelpoints(pc1(1),pc2(1),colNames_Day1_repRemoved(1))
hold on
for i = 2:2:8
    scatter(pc1(i:i+1),pc2(i:i+1),100,'filled')
    labelpoints(pc1(i:i+1),pc2(i:i+1),colNames_Day1_repRemoved(i:i+1))
    hold on
end

tree = linkage(day1_IBAQ_Scaled_RepRemoved');
figure;
dendrogram(tree,'labels',colNames_Day1_repRemoved)

day1_IBAQ_Scaled_RepRemoved(day1_IBAQ_Scaled_RepRemoved>=30)=NaN;
[r,pvalue] = corrplot(day1_IBAQ_Scaled_RepRemoved,'varNames',colNames_Day1_repRemoved,'rows','complete');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IBAQ Normalized Scalled to 100% - Day 3


[day3coeff,day3score,day3latent,day3tsquared,day3explained] = pca(day3_IBAQ_Scaled);
figure;
pc1 = day3coeff(:,1);
pc2 = day3coeff(:,2);
for i = 1:2:9
    scatter(pc1(i:i+1),pc2(i:i+1),100,'filled')
    labelpoints(pc1(i:i+1),pc2(i:i+1),colNames_Day3(i:i+1))
    hold on
end

tree = linkage(day3_IBAQ_Scaled');
figure;
dendrogram(tree,'labels',colNames_Day3)

day3_IBAQ_Scaled(day3_IBAQ_Scaled>=30)=NaN;
[r,pvalue] = corrplot(day3_IBAQ_Scaled,'varNames',colNames_Day3,'rows','complete');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('oldData_Day1.mat')
load('oldData_day3.mat')
load('oldData_day3_Rep34.mat')

%Original data from Robert - scaled to 100, no ibaq


[normcoeff,normscore,normlatent,normtsquared,day1OldExplained] = pca(oldData_Day1);
figure;
pc1 = normcoeff(:,1);
pc2 = normcoeff(:,2);
for i = 1:2:9
    scatter(pc1(i:i+1),pc2(i:i+1),100,'filled')
    labelpoints(pc1(i:i+1),pc2(i:i+1),colNames_Day1(i:i+1))
    hold on
end

tree = linkage(oldData_Day1');
figure;
dendrogram(tree,'labels',colNames_Day1)

oldData_Day1(oldData_Day1<=3)=NaN;
oldData_Day1(oldData_Day1>=30)=NaN;
[r,pvalue] = corrplot(oldData_Day1,'varNames',colNames_Day1,'rows','complete');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[normcoeff,normscore,normlatent,normtsquared,day3OldExplained] = pca(oldData_day3);
figure;
pc1 = normcoeff(:,1);
pc2 = normcoeff(:,2);
for i = 1:2:9
    scatter(pc1(i:i+1),pc2(i:i+1),100,'filled')
    labelpoints(pc1(i:i+1),pc2(i:i+1),colNames_Day3(i:i+1))
    hold on
end

tree = linkage(oldData_day3');
figure;
dendrogram(tree,'labels',colNames_Day3)

oldData_day3(oldData_day3<=3)=NaN;
oldData_day3(oldData_day3>=30)=NaN;
[r,pvalue] = corrplot(oldData_day3,'varNames',colNames_Day3,'rows','complete');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[normcoeff,normscore,normlatent,normtsquared,nromexplained] = pca(oldData_day3_Rep34);
figure;
pc1 = normcoeff(:,1);
pc2 = normcoeff(:,2);
for i = 1:2:9
    scatter(pc1(i:i+1),pc2(i:i+1),100,'filled')
    labelpoints(pc1(i:i+1),pc2(i:i+1),colNames_Day3(i:i+1))
    hold on
end

tree = linkage(oldData_day3_Rep34');
figure;
dendrogram(tree,'labels',colNames_Day3)

oldData_day3_Rep34(oldData_day3_Rep34<=3)=NaN;
oldData_day3_Rep34(oldData_day3_Rep34>=30)=NaN;
[r,pvalue] = corrplot(oldData_day3_Rep34,'varNames',colNames_Day3,'rows','complete');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fold change - day 1 and 3 together REMOVING BAD REPLICATE

load('day3_ibaq_scaled_matching.mat')
load('day1_ibaq_scaled_matching.mat')
load('colNames_Day1.mat')
load('colNames_Day3.mat')

% for i = 3:size(day1_ibaq_scaled_matching,2)
%     day1_fc(:,i-2) = day1_ibaq_scaled_matching(:,i)./day1_ibaq_scaled_matching(:,1);
% end
% 
% for i = 3:2:9
%     day3_fc(:,i-2) = day3_ibaq_scaled_matching(:,i)./day3_ibaq_scaled_matching(:,1);
%     day3_fc(:,i-1) = day3_ibaq_scaled_matching(:,i+1)./day3_ibaq_scaled_matching(:,2);
% end

for i = 3:10
    day1_fc(:,i-2) = day1_ibaq_scaled_matching(:,i)./day1_ibaq_scaled_matching(:,1);
end

for i = 3:10
    day3_fc(:,i-2) = day3_ibaq_scaled_matching(:,i)./mean(day3_ibaq_scaled_matching(:,1:2),2);
end

day1_log2fc = log2(day1_fc);
day3_log2fc = log2(day3_fc); 


day1_3_log2fc = [day1_log2fc day3_log2fc];
day1_3_names = [colNames_Day1(3:end) colNames_Day3(3:end)];

%%%%%%%%%%%%%%%%%%%%%%


[fccoeff,fcscore,fclatent,fctsquared,fcexplained] = pca(day1_3_log2fc);
figure;
pc1 = fccoeff(:,1);
pc2 = fccoeff(:,2);
for i = 1:2:15
    scatter(pc1(i:i+1),pc2(i:i+1),100,'filled')
    labelpoints(pc1(i:i+1),pc2(i:i+1),day1_3_names(i:i+1))
    hold on
end

tree = linkage(day1_3_log2fc');
figure;
dendrogram(tree,'labels',day1_3_names)

[r,pvalue] = corrplot(day1_3_log2fc,'varNames',day1_3_names,'rows','complete');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Fold change - day 1 and 3 together keep all replicates

load('day3_ibaq_scaled_matching.mat')
load('day1_ibaq_scaled_matching.mat')
load('colNames_Day1.mat')
load('colNames_Day3.mat')

% for i = 3:2:9
%     day1_fc(:,i-2) = day1_ibaq_scaled_matching(:,i)./day1_ibaq_scaled_matching(:,1);
%     day1_fc(:,i-1) = day1_ibaq_scaled_matching(:,i+1)./day1_ibaq_scaled_matching(:,2);
% end
% 
% for i = 3:2:9
%     day3_fc(:,i-2) = day3_ibaq_scaled_matching(:,i)./day3_ibaq_scaled_matching(:,1);
%     day3_fc(:,i-1) = day3_ibaq_scaled_matching(:,i+1)./day3_ibaq_scaled_matching(:,2);
% end


for i = 3:10
    day1_fc(:,i-2) = day1_ibaq_scaled_matching(:,i)./mean(day1_ibaq_scaled_matching(:,1:2),2);
end

for i = 3:10
    day3_fc(:,i-2) = day3_ibaq_scaled_matching(:,i)./mean(day3_ibaq_scaled_matching(:,1:2),2);
end



day1_log2fc = log2(day1_fc);
day3_log2fc = log2(day3_fc); 


day1_3_log2fc = [day1_log2fc day3_log2fc];
day1_3_names = [colNames_Day1(3:end) colNames_Day3(3:end)];

%%%%%%%%%%%%%%%%%%%%%%

[r,pvalue] = corrplot(day1_3_log2fc,'varNames',day1_3_names,'rows','complete');

[fccoeff,fcscore,fclatent,fctsquared,fcexplained] = pca(day1_3_log2fc);
figure;
pc1 = fccoeff(:,1);
pc2 = fccoeff(:,2);
for i = 1:2:15
    scatter(pc1(i:i+1),pc2(i:i+1),100,'filled')
    labelpoints(pc1(i:i+1),pc2(i:i+1),day1_3_names(i:i+1))
    hold on
end

tree = linkage(day1_3_log2fc');
figure;
dendrogram(tree,'labels',day1_3_names)
