function [idx,elapsed_time, Wnbc2, basins] = nbc2(X,fhandle,r,b)
% Reference:
%
% Mike Preuss.
% Improved Topological Niching for Real-Valued Global Optimization,
% EvoApplications 2012, LNCS 7248, pp. 386--395, 2012.


[n,d] = size(X);

if nargin < 4, 
	b = (-4.69e-4*d^2  + 0.0263*d + 3.66/d - 0.457)*log10(n) + ...
		7.51e-4*d^2 + 0.0421*d - 2.26/d + 1.83;
end

if nargin < 3, r = 2; end;


global initial_flag;
initial_flag = 0;

fx = zeros(n,1);
for i=1:n,
	fx(i) = fhandle(X(i,:));
end
assert(sum(isnan(fx))==0, 'NaN in fx');
assert(sum(isinf(fx))==0, 'Infinite in fx');

start_time = tic;
% distance matrix
W = spalloc(n,n,2*n);
outgoing = zeros(n,2);
for i=1:n,
	index = find(fx > fx(i));
	if length(index) > 0,
		[ds,j] = min(pdist2( X(i,:), X(index,:) ));

		% To avoid problems arising from identical observations
		% which have ds=0 and therefore appear as non-entries in W
		ds = max(ds, eps);

		% Columns contain incoming connections - ONLY -
		W(i, index(j)) = ds;

		% node to which i is connected
		outgoing(i,1) = index(j);
		% length of outgoing connection of node i
		outgoing(i,2) = ds;
	end
end

% Rule 1
% Delete all edges that are longer than r times the average length
av = mean(nonzeros(W));
[i,j] = find(W > r*av);
W(i,j) = 0;

if d > 2,
	% Rule 2: Delete outgoing connections 
	for i=1:n,
		% Incoming edges
		inc = nonzeros(W(:,i));
		% For points with at least 3 incoming and 1 outgoing edge 
		% (second part of rule is irrelevant since I am deleting outgoing edge)
		if length(inc) > 2 && outgoing(i,2)/median(inc) > b,
			W(i, outgoing(i,1)) = 0;
		end
	end
end

% Find connected components in graph
Wnbc2 = adjacency( graph(max(W,W')) );
basins = conncomp(graph(max(W,W')), 'OutputForm', 'cell');

idx = zeros(n,1);
for i=1:length(basins),
	idx( basins{i} ) = i;
end

elapsed_time = toc(start_time);
