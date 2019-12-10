function [W,M] = adjacency(A, X, y, fhandle, budget)
%Estimate affinity matrix
%function [W,M] = affinity(A, X, y, fhandle, budget)
%
%Inputs:
% A: {0,1} Affinity matrix defining edges in graph
% X: Matrix storing points row-wise
% y: Function values at these points
% fhandle: Function handle to objective function
% B: Budget of function evaluations
%
%Outputs:
% W: Weigthed affinity matrix
% M: Matrix indicating local minima on edge weights: 
%	M(i,j) = 0; edge does not exist
%	M(i,j) = 1; local minimum found
%	M(i,j) = 2; no function evaluations on link
%	M(i,j) = 3; no local minimum detected

assert(size(A,1)==size(X,1), 'Data matrix and Adjacency matrix are inconsistent');

if nargin<5,
	error('ADJACENCY: All inputs needs to be specified\n');
end

n = size(A,1);
W = spalloc(n,n, n+nnz(A));
M = W;

if ~isempty(y)
	% Set the diagonal of similarity matrix in vector format
	W(1:(n+1):(n*n)) = y';
end


% Estimate edge lengths and sum of lengths of all edges
D = spalloc(n,n, nnz(A));
total_dist = 0;
for i=1:n,
	% find observations with which i is connected
	links = find(A(i,i+1:n)>0) + i;
	% Compute Euclidean distances
	dst = pdist2(X(i,:), X(links,:));

	D(i, links) = dst;
	total_dist = total_dist + sum(dst);
end
% function evaluations per distance unit
fpd = budget/total_dist;

% Identify constant b^\star in the paper (see Eq.(9))
a = 1;
s = sum(round(D(:) * fpd*a));
% bound the range of solution
if s < budget,
	b = 2;
	while sum(round(D(:) * fpd*b)) < budget,
		b = 2*b;
	end
elseif s > budget,
	b = 0.5; 
	while sum(round(D * fpd * b)) > budget,
		b = 0.5*b;
	end
	a = b;
	b = 1;
end
% required for the case when s==budget in first iteration
c = a;

% Bisection
while abs(s-budget) > 0 & abs(b-a)> sqrt(eps),
	c = 0.5*(a + b);
	s = sum(round(D(:) * fpd * c));
	if s > budget,
		b = c;
	else
		a = c;
	end
end
%%%%%%%%%%% Number of fcn evaluations per edge
D = round(D*fpd*c);

% DEBUGGING
%if sum(D(:)) ~= budget,
%	fprintf('sn_similarity_budget: difference in actual versus specified budget: %i ', sum(D(:)) - budget);
%	%keyboard;
%end


for i = 1:n-1,
	% Evaluate similarity to connected nodes only
	for j = find(A(i,i+1:n)>0) + i,
		fmin = min(y(i), y(j));
		% Allocate budget on specific edge
		a = linspace(0,1, D(i,j) + 2);

		if length(a) == 2,
			% no function evaluations on this link
			M(i,j) = 2;
		else
			% no local minimum found during search
			M(i,j) = 3;
		end
		
		% initialise
		f = zeros(length(a),1);
		f(1) = y(i);
		f(end) = y(j);
		for l=2:length(a)-1,
			f(l) = fhandle( (1-a(l))*X(i,:) + a(l)* X(j,:) );
			if f(l) < fmin & abs(f(l)-fmin)>sqrt(eps),
				fmin = f(l);
				% local minimum found
				M(i,j)=1;
			end
		end

		W(i,j) = fmin;
		W(j,i) = W(i,j);

		M(j,i) = M(i,j);
	end
end
