function Wsym = slink(X, Wsym, comp),
%Single Linkage on Connectivity Graph
%function Wsym = slink(X, Wsym, comp),
%
%Inputs:
% X: data matrix containing points (vertices)
% Wsym: affinity matrix describing graph
% comp: cell array each element of which is a row vector 
%	of the indices (row numbers) of observations that belong
%	to the given connected component of the graph 
%
%Output:
% W: Affinity matrix

if length(comp)==1,
	return;
end

K = length(comp);
% single linkage distance between connected components
Z = inf*ones(K,K);
% nodes in connected components achieving min distance
L = zeros(K,K);
for i=1:K-1,
	for j=(i+1):K,
		M = pdist2(X(comp{i},:),  X(comp{j},:));

		[val, k] = min(M(:));
		Z(i,j) = val;
		Z(j,i) = val;

		% Location of minimiser in M
		[pi, pj] = ind2sub(size(M), k);
		% id (row in data matrix X) of node in comp{i} closest to comp{j}
		L(i,j) = comp{i}(pi);
		% id (row in data matrix X) of node in comp{j} closest to comp{i}
		L(j,i) = comp{j}(pj);
	end
end

% Second pass: Iteratively connect closest (in terms of single linkage) components
for r = 1:K-1,
	% minimum distance between components
	[val, k] = min(Z(:));
	% components that minimise this distance
	[i, j] = ind2sub(size(Z),k);
	% this helps with the overwriting later on
	if i > j,
		cp = i;
		i = j;
		j = cp;
	end

	% add edge in adjacency matrix
	Wsym(L(i,j), L(j,i)) = 1;
	Wsym(L(j,i), L(i,j)) = 1;

	% identify minimum distance of merged component to all others
	[minD, loc] = min(Z(:,[i j]),[],2);

	% Update Z: 
	Z(:,i) = minD;
	Z(i,:) = minD;
	%	strictly speaking Z(i,j) need not be set
	Z([i j],i) = inf;
	%	cancel out j-th column	
	Z(:,j) = inf;
	Z(j,:) = inf;

	% Update nearest neighbor information
	% Matrix M(i,j) = 1 if observation from comp{i} is closest to comp{j}
	M = sparse(loc,1:K,ones(K,1));
	L(i,:) = max( [L(i,:); L(j,:)].*M );
	L(i,i) = 0;
	L(j,:) = 0;
	L(:,j) = 0;

	% merge i and j and then empty j
	comp{i} = [comp{i},comp{j}];
	comp{j} = [];
end
