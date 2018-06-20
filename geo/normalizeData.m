function [normData] = normalizeData(tkiData,doxData)


jointData = [tkiData doxData]; 

%Scale data to mean 0, std 1 
for row = 1:size(jointData,1)
    centeredData(row,:) = jointData(row,:)-nanmean(jointData(row,:));
end

for row = 1:size(jointData,1)
    normData(row,:) = centeredData(row,:)./nanstd(jointData(row,:));
end

normData(normData==inf)=NaN;
normData(normData==-inf)=NaN;
normData(isnan(normData))=0;

for i = 1:size(normData,1)
    if ~isnan(normData(i,end))
        normData(i,:) = normData(i,:);
    end
end
% normData = normData2;
