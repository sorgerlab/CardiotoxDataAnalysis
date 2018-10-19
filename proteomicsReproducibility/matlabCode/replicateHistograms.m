set(0,'DefaultFigureVisible','off');

load('2016_day3_allReplicates.mat')
% load('day3_allReplicates.mat')
load('replicate_labels.mat')


%1 v 2 (Technical 1)
%3 v 4 (Technical 2)
%1 v 3 (Biological 1)
%2 v 4 (Biological 2)
%Potentially 1 v 4, 2 v 3

%Need to remove 0 rows
day3_allReplicates(day3_allReplicates <=1) = NaN;

%Technical Replicates
for i = 1:2:size(day3_allReplicates,2)
    relicateFoldChange = log2(day3_allReplicates(:,i)./day3_allReplicates(:,i+1));
%     subplot(10,1,1)
    figure
    histogram(relicateFoldChange)
    xlim([-1 1]);
    title([replicate_labels{i} 'vs' replicate_labels{i+1}],'Interpreter', 'none');
    ylabel('Frequence','fontsize',12,'fontweight','b','fontname','Arial');
    xlabel('Log2(Rep 1/Rep 2)','fontsize',12,'fontweight','b','fontname','Arial');
    set(gca,'fontsize',10,'fontname','Arial');
    set(gcf,'color','w');
    meanFC = sprintf('%.3f',nanmean(relicateFoldChange));
    varFC =sprintf('%.3f',nanstd(relicateFoldChange)^2);
    annotation('textbox',[.65 .7 .1 .2],'String',{['Mean:' meanFC] , ['Var:' varFC]},'EdgeColor','none','fontsize',16,'fontname','Arial')
%     filename = ['Figures/2017/Histograms/TechnicalReplicates/' [replicate_labels{i} '_' replicate_labels{i+1}] '.png'];
    filename = ['Figures/2016/Histograms/TechnicalReplicates/' [replicate_labels{i} '_' replicate_labels{i+1}] '.png'];
    print(filename, '-dpng','-r600')
   
end


%Subset of Biological Replicates
for i = [1,2,5,6,9,10,13,14,17,18]
    relicateFoldChange = log2(day3_allReplicates(:,i)./day3_allReplicates(:,i+2));
%     subplot(10,1,1)
    figure
    histogram(relicateFoldChange)
    xlim([-1 1]);
    title([replicate_labels{i} 'vs' replicate_labels{i+2}],'Interpreter', 'none');
    ylabel('Frequence','fontsize',12,'fontweight','b','fontname','Arial');
    xlabel('Log2(Rep 1/Rep 2)','fontsize',12,'fontweight','b','fontname','Arial');
    set(gca,'fontsize',10,'fontname','Arial');
    set(gcf,'color','w');
    meanFC = sprintf('%.3f',nanmean(relicateFoldChange));
    varFC =sprintf('%.3f',nanstd(relicateFoldChange)^2);
    annotation('textbox',[.65 .7 .1 .2],'String',{['Mean:' meanFC] , ['Var:' varFC]},'EdgeColor','none','fontsize',16,'fontname','Arial')
%     filename = ['Figures/2017/Histograms/BiologicalReplicates/' [replicate_labels{i} '_' replicate_labels{i+2}] '.png'];
    filename = ['Figures/2016/Histograms/BiologicalReplicates/' [replicate_labels{i} '_' replicate_labels{i+2}] '.png'];
    print(filename, '-dpng','-r600')

end 

%Negative control - dmso vs sorafenib
%Subset of Biological Replicates
for i = 1:4
    relicateFoldChange = log2(day3_allReplicates(:,i)./day3_allReplicates(:,i+12));
%     subplot(10,1,1)
    figure
    histogram(relicateFoldChange)
    xlim([-1 1]);
    title([replicate_labels{i} 'vs' replicate_labels{i+12}],'Interpreter', 'none');
    ylabel('Frequence','fontsize',12,'fontweight','b','fontname','Arial');
    xlabel('Log2(Rep 1/Rep 2)','fontsize',12,'fontweight','b','fontname','Arial');
    set(gca,'fontsize',10,'fontname','Arial');
    set(gcf,'color','w');
    meanFC = sprintf('%.3f',nanmean(relicateFoldChange));
    varFC =sprintf('%.3f',nanstd(relicateFoldChange)^2);
    annotation('textbox',[.65 .7 .1 .2],'String',{['Mean:' meanFC] , ['Var:' varFC]},'EdgeColor','none','fontsize',16,'fontname','Arial')
%     filename = ['Figures/2017/Histograms/NegativeControl/' [replicate_labels{i} '_' replicate_labels{i+12}] '.png'];
    filename = ['Figures/2016/Histograms/NegativeControl/' [replicate_labels{i} '_' replicate_labels{i+12}] '.png'];
    print(filename, '-dpng','-r600')

end 
