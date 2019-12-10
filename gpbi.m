function [idx, noise, fc, ctr, S, elapsed] = gpbi(X, f, k, B, C)
%function [idx, elapsed] = gpbi(X, f, k, B, C)
%
%Inputs:
% X: Matrix whose rows contain points in the domain of f
% f: Function handle to objective function @(x)(f(x)) where x is a row vector
% k: Number of k-Nearest Neighbours used to create graph
% B: Budget of function evaluations to estimate graph affinity matrix
% C: Number of function evaluations for each merge operation during pruning

if ~isdeployed,
	addpath(genpath('src'));
end

n = size(X,1);
y = zeros(n,1);
for i=1:n,
	y(i) = f( X(i,:) );
end

startC = tic;
% Connectivity Graph (connect to k-Nearest Neighbours and Nearest Better Neighbour)
W = cgraph(X,y,k);
comp = conncomp(graph(W,'upper'), 'OutputForm', 'cell');

% if graph is disconnected perform single linkage
if length(comp) > 1,
	fprintf('Connected components before SL %i\n', length(comp));
	W = slink(X,W,comp);
end

% Estimate edge weights
[Wsn,Msn] = adjacency(W, X, y, f, B);

% Graph Partitioning, Pruning and Noise allocation
[idx,noise,fc,ctr,S] = gpart(Wsn, Msn, X, f, C);
elapsed = toc(startC);
