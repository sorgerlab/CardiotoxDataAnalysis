load('Day3R1.mat') %D3_R1~rq_131_sn_scaled
load('Day3R2.mat') %D3_R2~rq_128c_sn_scaled
load('colNames.mat')

combinedData = [Day3R1 Day3R2];    

%Need to remove 0 rows
Day3R1(Day3R1 <=1) = NaN;
Day3R2(Day3R1 <=1) = NaN;
Day3R1(Day3R2 <=1) = NaN;
Day3R2(Day3R2 <=1) = NaN;

%Reorganize data, so there is one array for each drug
%Each array contains two biological replicates, each of which contain two
%technical replicates
%Take mean of Technical Replicate 1 from Biological Replicate 1 and 2
%Take mean of Technical Replicate 2 from Biological Replicate 1 and 2
%Treat as final two replicates for comparison 
dmsoData = [nanmean([Day3R1(:,1), Day3R2(:,1)],2), nanmean([Day3R1(:,2), Day3R2(:,2)],2)];
erlData = [nanmean([Day3R1(:,3), Day3R2(:,3)],2), nanmean([Day3R1(:,4), Day3R2(:,4)],2)];
lapData = [nanmean([Day3R1(:,5), Day3R2(:,5)],2), nanmean([Day3R1(:,6), Day3R2(:,6)],2)];
sorData = [nanmean([Day3R1(:,7), Day3R2(:,7)],2), nanmean([Day3R1(:,8), Day3R2(:,8)],2)];
sunData = [nanmean([Day3R1(:,9), Day3R2(:,9)],2), nanmean([Day3R1(:,10), Day3R2(:,10)],2)];

%Take log2 fold change of data
dmso = log2(dmsoData(:,1)./dmsoData(:,2));
erl = log2(erlData(:,1)./erlData(:,2));
lap = log2(lapData(:,1)./lapData(:,2));
sor = log2(sorData(:,1)./sorData(:,2));
sun = log2(sunData(:,1)./sunData(:,2));

%Plot log2 FC for replicates of each drug
figure
subplot(5,1,1)
histogram(dmso)
xlim([-1 1]);
title('DMSO');
ylabel('Frequence','fontsize',12,'fontweight','b','fontname','Arial');
xlabel('Log2(Rep 1/Rep 2)','fontsize',12,'fontweight','b','fontname','Arial');
set(gca,'fontsize',10,'fontname','Arial');
set(gcf,'color','w');

subplot(5,1,2)
histogram(erl)
xlim([-1 1]);
title('Erlotinib');
ylabel('Frequence','fontsize',12,'fontweight','b','fontname','Arial');
xlabel('Log2(Rep 1/Rep 2)','fontsize',12,'fontweight','b','fontname','Arial');
set(gca,'fontsize',10,'fontname','Arial');
set(gcf,'color','w');

subplot(5,1,3)
histogram(lap)
xlim([-1 1]);
title('Lapatinib');
ylabel('Frequence','fontsize',12,'fontweight','b','fontname','Arial');
xlabel('Log2(Rep 1/Rep 2)','fontsize',12,'fontweight','b','fontname','Arial');
set(gca,'fontsize',10,'fontname','Arial');
set(gcf,'color','w');

subplot(5,1,4)
histogram(sor)
xlim([-1 1]);
title('Sorafenib');
ylabel('Frequence','fontsize',12,'fontweight','b','fontname','Arial');
xlabel('Log2(Rep 1/Rep 2)','fontsize',12,'fontweight','b','fontname','Arial');
set(gca,'fontsize',10,'fontname','Arial');
set(gcf,'color','w');

subplot(5,1,5)
histogram(sun)
xlim([-1 1]);
title('Sunitinib');
ylabel('Frequence','fontsize',12,'fontweight','b','fontname','Arial');
xlabel('Log2(Rep 1/Rep 2)','fontsize',12,'fontweight','b','fontname','Arial');
set(gca,'fontsize',10,'fontname','Arial');
set(gcf,'color','w');