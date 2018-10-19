function [normData] = normalizeData(rawData)

% Scale data to mean 0, std 1 
centeredData = zeros(size(rawData));

for col = 1:size(rawData,2)
    centeredData(:,col) = rawData(:,col)-nanmean(rawData(:,col));
end

normData = zeros(size(rawData));
for col = 1:size(rawData,2)
    normData(:,col) = centeredData(:,col)./nanstd(centeredData(:,col));
end

normData(normData==inf)=NaN;
normData(normData==-inf)=NaN;
normData(any(isnan(normData), 2), :) = [];
