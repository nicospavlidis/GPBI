% Script that illustrates the assignment of solutions produced by NEA2 algorithm
% for CEC2013 MMO competition

clear all;

% Specify CEC function number
cecF = 7;
K = 4;
% fcn evaluation proportion
b = 0.5;
C = 10;

% Reset flag (necessary for Composition Functions)
global initial_flag;
initial_flag = 0;
fh = @(x)(niching_func(x,cecF));

addpath('compete');
addpath(genpath('cec13mmo'));

X = load(sprintf('nea2solutions/nea2_%i_X_unique.txt',cecF));
labels = load(sprintf('nea2solutions/label_%i.txt', cecF));
fprintf('Data loaded successfully\n');

[n,d] = size(X);

T = table([],[],[],[],[],[],[],[],'VariableNames',{'fcn','method', 'K', 'budget', 'nclust','nmi', 'ari','purity'});

% NBC2
[idx,elapsed_time, Wnbc2] = nbc2(X,fh);
pn = performance(idx, labels);
T = [T; {cecF, 'nbc2',0, 0, length(unique(idx)), pn.NMI, pn.AdjRand, pn.Purity}];
fprintf('NBC2 complete\n');

% BAPD
[idx,elapsed_time] = bapd(X,fh);
pn = performance(idx, labels);
T = [T; {cecF, 'bapd',0, 0, length(unique(idx)), pn.NMI, pn.AdjRand, pn.Purity}];
fprintf('BAPD complete\n');

% SEEDS
[idx,elapsed_time] = seeding(X,fh,[],cecF);
pn = performance(idx, labels);
T = [T; {cecF,'seed',0, 0, length(unique(idx)), pn.NMI, pn.AdjRand, pn.Purity}];
fprintf('SEEDS complete\n');

% GPBI
B = floor(b*K*n);
[idx,noise,~,~,S,myclock] = gpbi(X, fh, K, B, C);
pn = performance(idx, labels);
T = [T; {cecF,'gpbi',K, B, max(idx), pn.NMI, pn.AdjRand, pn.Purity}];
fprintf('GPBI complete\n');

writetable(T,sprintf('perf_cecF%i.txt',cecF));
