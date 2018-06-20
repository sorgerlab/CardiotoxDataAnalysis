% To set axes:
set(gca,'XLim',[1 5],'XTick',1.5:4.5,'XTickLabel',xticklabels,'Fontname','Arial','FontWeight','Bold','Fontsize',7)
XTick can take a custom array of numbers, instead of being evenly spaced between upper and lower limits, such as [1,3,8]
XTickLabel can take an array of numbers or cell array of strings, such as xticklabels={'6 h','24 h','72 h','168 h'};

% To set figure size, when first creating a figure:
figure('Units','centimeters', 'Position', [100, 100, 15.325, 13.031]);

% To set color bar position:
% In this case I found the location of the axes for the 20th subplot I created, so the 5th row, 4th column of the 7x4 crossplot array
% This was the closest I could get to the desired positioning automatically
% From there I played with the numbers to get the exact size and position I wanted

hp4 = get(subplot(7,4,20),'Position');
h=colorbar('Position', [hp4(1)+hp4(3)+0.015  hp4(2)-.175  0.0275  hp4(2)+hp4(3)*2.1]);
set(h,'Fontname','Arial','FontWeight','Bold','Fontsize',7);
    
% Also potentially useful, to print to a file of a certain name, file type, and resolution:    
print('filename.svg', '-Painters', '-dsvg','-r600')