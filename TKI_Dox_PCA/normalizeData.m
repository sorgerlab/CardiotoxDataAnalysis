function [normData] = normalizeData(tkiData,doxData)

jointData = [tkiData doxData]; 

% Scale data to mean 0, std 1 
centeredData = zeros(size(jointData));
for row = 1:size(jointData,1)
    centeredData(row,:) = jointData(row,:)-nanmean(jointData(row,:));
end

normData = zeros(size(jointData));
for row = 1:size(jointData,1)
    normData(row,:) = centeredData(row,:)./nanstd(jointData(row,:));
%     mean(normData(row,:))
%     std(normData(row,:))
end

normData(normData==inf)=NaN;
normData(normData==-inf)=NaN;
normData(any(isnan(normData), 2), :) = [];
