set(0,'DefaultFigureVisible','on');

load('synapse_day3_allReplicates.mat')
% load('day3_allReplicates.mat')
load('replicate_labels.mat')

logData = log10(day3_allReplicates);
logData(logData == inf) = NaN;
logData(logData == -inf) = NaN;


j=1
filenames = {'dmso_correlations.png','erl_correlations.png','lap_correlations.png','soraf_correlations.png','sunit_correlations.png'};
for i = 1:4:size(logData,2)
    
    figure;
    [r2,pvalue] = corrplot(logData(:,i:i+3),'varNames',replicate_labels(i:i+3));

    hfig = gcf;
    haxes = findobj(hfig, 'Type', 'Axes');
    set(haxes,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
    set(0,'defaultfigurecolor',[1 1 1])

%     print(['Figures/2017/Correlations/' filenames{j}], '-Painters', '-dpng','-r600')
    print(['Figures/2016/Correlations/' filenames{j}], '-Painters', '-dpng','-r600')
    j=j+1;
end













