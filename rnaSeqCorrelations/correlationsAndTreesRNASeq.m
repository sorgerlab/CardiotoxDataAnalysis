load('allDmso.mat')
load('allErl.mat')
load('allLapat.mat')
load('allSoraf.mat')
load('allSunit.mat')
load('dmsoNames.mat')
load('erlNames.mat')
load('lapatNames.mat')
load('sorafNames.mat')
load('sunitNames.mat')

replicates = {'Rep 1', 'Rep 2', 'Rep 3'};

% Pearson correlation of TMM-normalized counts data between replicates of DMSO treatment at 24 hour and Sorafenib treatment at 10µM and 24 hours were plotted. 

%%%%Correlation%%%%%
%Current Names
% {'S1C1','S1C2','S1C3','S2B1','S2B2','S2B3','S2C1','S2C2','S2C3','S2D1','S2D2','S2D3','S3C1','S3C2','S3C3','S4C1','S4C3'}
%Want 10 µM and 24 hours, so high dose, mid time
%So So 2 - D replicates - 'S2D1','S2D2','S2D3'
% dmso 4:6
%Columns 10, 11, 12 for soraf.
% Sunit 10:11
% Erl 10:12
% lapat 10:12



dmso_10u_24hr = dmsoRNA(:,4:6);
logdmso = log10(dmso_10u_24hr);
logdmso(logdmso == inf) = NaN;
logdmso(logdmso == -inf) = NaN;

figure;
[r2,pvalue] = corrplot(logdmso,'varNames',namesToUse);

hfig = gcf;
haxes = findobj(hfig, 'Type', 'Axes');
set(haxes,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
set(0,'defaultfigurecolor',[1 1 1])

fn1 = 'dmso.svg';
print(fn1, '-Painters', '-dsvg','-r600')




erl_10u_24hr = erlotinibRNA(:,10:12);
logerl = log10(erl_10u_24hr);
logerl(logerl == inf) = NaN;
logerl(logerl == -inf) = NaN;

[r2,pvalue] = corrplot(logerl,'varNames',namesToUse);

hfig = gcf;
haxes = findobj(hfig, 'Type', 'Axes');
set(haxes,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
set(0,'defaultfigurecolor',[1 1 1])

fn1 = 'erl.svg';
print(fn1, '-Painters', '-dsvg','-r600')



lapat_10u_24hr = lapatinibRNA(:,10:12);
loglapat = log10(lapat_10u_24hr);
loglapat(loglapat == inf) = NaN;
loglapat(loglapat == -inf) = NaN;

[r2,pvalue] = corrplot(loglapat,'varNames',namesToUse);


hfig = gcf;
haxes = findobj(hfig, 'Type', 'Axes');
set(haxes,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
set(0,'defaultfigurecolor',[1 1 1])

fn1 = 'lap.svg';
print(fn1, '-Painters', '-dsvg','-r600')


% figure;
soraf_10u_24hre = sorafenibRNA(:,10:12);
logSoraf = log10(soraf_10u_24hre);
logSoraf(logSoraf == inf) = NaN;
logSoraf(logSoraf == -inf) = NaN;

 [r2,pvalue] = corrplot(logSoraf,'varNames',namesToUse);

 hfig = gcf;
haxes = findobj(hfig, 'Type', 'Axes');
set(haxes,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
set(0,'defaultfigurecolor',[1 1 1])

fn1 = 'sor.svg';
print(fn1, '-Painters', '-dsvg','-r600')



sunit_10u_24hr = sunitinibRNA(:,10:11);
logsunit = log10(sunit_10u_24hr);
logsunit(logsunit == inf) = NaN;
logsunit(logsunit == -inf) = NaN;

[r2,pvalue] = corrplot(logsunit,'varNames',namesToUse([1 3]));

hfig = gcf;
haxes = findobj(hfig, 'Type', 'Axes');
set(haxes,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
set(0,'defaultfigurecolor',[1 1 1])

fn1 = 'sun.svg';
print(fn1, '-Painters', '-dsvg','-r600')






%%%%%%%%%%%%%%%%
%Hierarchical clusters
% 
% tree = linkage(allDmso');
% figure;
% [h,nodes]=dendrogram(tree,'labels',dmsoNames)
% 
% tree = linkage(allSoraf');
% figure;
% [h,nodes]=dendrogram(tree,'labels',sorafNames)
% 
% tree = linkage(allLapat');
% figure;
% [h,nodes]=dendrogram(tree,'labels',lapatNames)
