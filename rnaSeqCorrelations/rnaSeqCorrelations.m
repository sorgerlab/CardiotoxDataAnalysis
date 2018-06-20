function [r,pvalue] = rnaSeqCorrelations()

load('dmsoRNA.mat')
load('erlotinibRNA.mat')
load('lapatinibRNA.mat')
load('sorafenibRNA.mat')
load('sunitinibRNA.mat')


replicates = {'Rep 1', 'Rep 2', 'Rep 3'};

%Plot earson correlation of TMM-normalized counts data between replicates 
%Using 10µM and 24 hours treatment conditions
%Corresponding Data in columns 4:6 in DMSO data, 10:12 in TKI Data

%%%%% DMSO %%%%%

%Extract and log transform counts data
dmso_10u_24hr = dmsoRNA(:,4:6);
logdmso = log10(dmso_10u_24hr);
logdmso(logdmso == inf) = NaN;
logdmso(logdmso == -inf) = NaN;

%Plot correlations between replicates
figure;
[r,pvalue] = corrplot(logdmso,'varNames',replicates);

%Set axis label options
hfig = gcf;
haxes = findobj(hfig, 'Type', 'Axes');
set(haxes,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
set(0,'defaultfigurecolor',[1 1 1])

%%%%% Erlotinib %%%%%

%Extract and log transform counts data
erl_10u_24hr = erlotinibRNA(:,10:12);
logerl = log10(erl_10u_24hr);
logerl(logerl == inf) = NaN;
logerl(logerl == -inf) = NaN;

%Plot correlations between replicates
figure;
[r,pvalue] = corrplot(logerl,'varNames',replicates);

%Set axis label options
hfig = gcf;
haxes = findobj(hfig, 'Type', 'Axes');
set(haxes,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
set(0,'defaultfigurecolor',[1 1 1])

%%%%% Lapatinib %%%%%

%Extract and log transform counts data
lapat_10u_24hr = lapatinibRNA(:,10:12);
loglapat = log10(lapat_10u_24hr);
loglapat(loglapat == inf) = NaN;
loglapat(loglapat == -inf) = NaN;

%Plot correlations between replicates
figure;
[r,pvalue] = corrplot(loglapat,'varNames',replicates);

%Set axis label options
hfig = gcf;
haxes = findobj(hfig, 'Type', 'Axes');
set(haxes,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
set(0,'defaultfigurecolor',[1 1 1])

%%%%% Sorafenib %%%%%

%Extract and log transform counts data
soraf_10u_24hre = sorafenibRNA(:,10:12);
logSoraf = log10(soraf_10u_24hre);
logSoraf(logSoraf == inf) = NaN;
logSoraf(logSoraf == -inf) = NaN;

%Plot correlations between replicates
figure;
[r,pvalue] = corrplot(logSoraf,'varNames',replicates);

%Set axis label options
hfig = gcf;
haxes = findobj(hfig, 'Type', 'Axes');
set(haxes,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
set(0,'defaultfigurecolor',[1 1 1])

%%%%% Sunitinib %%%%%

%Extract and log transform counts data
sunit_10u_24hr = sunitinibRNA(:,10:11);
logsunit = log10(sunit_10u_24hr);
logsunit(logsunit == inf) = NaN;
logsunit(logsunit == -inf) = NaN;

%Plot correlations between replicates
figure;
[r,pvalue] = corrplot(logsunit,'varNames',replicates([1 3]));

%Set axis label options
hfig = gcf;
haxes = findobj(hfig, 'Type', 'Axes');
set(haxes,'Fontname','Arial','FontWeight','Bold','Fontsize',12)
set(0,'defaultfigurecolor',[1 1 1])





