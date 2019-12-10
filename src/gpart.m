function [idx, noise, fcount, ctr, A] = gpart(Wsn,Msn,X,fhandle, pB)
%Graph Partitioning, Pruning and Noise Point Allocation
%function [idx, noise, fcount, ctr, A] = gpart(W,M,X,f)
%
%Inputs:
% W: Weighted adjacency matrix 
% M: Matrix indicating presence of minima
% X: Matrix storing points in rows
% f: Function handle to objective function
% pB: Function evaluation budget for Pruning

if ~isdeployed,
	addpath(genpath('src'));
end


A = Wsn;
% use (full) to avoid zeros appearing as all zero sparse matrices
fx = full(diag(A));
n = size(A,1);

% Maximum spanning tree (A) = - Minimum Spaning Tree(-A)
% graph which consists of 2 connected components 
A = -mst(-A);

[maxf, mode] = max(fx);
% initial split value is max(f(X(index,:))), index := points in current cluster
nd = struct('index',[1:n]', 'level', min(A(:)), 'split', maxf, ...
		'mode', mode,'depth', 0, 'link',[0,0]);
ctr = tree(nd);

l = 1;
while (1),
	% breadth first cluster splitting (since every node added to ctr is
	% assigned to the back of the Node{} cell array)
	leaf = ctr.findleaves();
	while (l <= length(leaf)) & (length( ctr.Node{leaf(l)}.index ) == 1),
		l = l+1;
	end
	if l > length(leaf),
		break;
	end
	leaf = leaf(l);

	% Indices of observations assigned to current leaf
	index = ctr.Node{leaf}.index;

	% Upper triangular of Adjacency matrix excluding diagonal
	U = triu(A(index, index), 1);

	% DEBUGGING
	if nnz(U)==0,
		error('Empty subgraph considered. This should never have happened\n');
		%keyboard;
	end

	% iterate over all links in connected component
	for i = 1:nnz(U),
		% identify smallest weight
		mW = min(nonzeros(U));
		[ind_i, ind_j] = find(U == mW);
		
		% in case minimum weight is not unique remove smallest subgraph
		if length(ind_i) > 1, 
			% Remove edge that results in smallest minimum size of C{1},C{2}
			minSize = inf;
			jstar = 0;
			for j=1:length(ind_i),
				% Remove connection
				U(ind_i(j), ind_j(j)) = 0;

				% link gives id (row numbers) of observations connected by
				% arch that is removed last
				%link = [index(ind_i(j)), index(ind_j(j))];

				% there can be only two connected components after removing one link
				C = conncomp(graph(U,'upper'), 'OutputForm', 'cell');
				mS = min(length(C{1}), length(C{2}));
				if mS < minSize,
					minSize = mS;
					jstar = j;
				end
				U(ind_i(j), ind_j(j)) = mW;
			end
			ind_i = ind_i(jstar);
			ind_j = ind_j(jstar);
		end

		% Remove edge that results in smallest minimum size of C{1},C{2}
		U(ind_i(1), ind_j(1)) = 0;

		% link gives id (row numbers) of observations connected by
		% edge that is removed last
		link = [index(ind_i(1)), index(ind_j(1))];

		C = conncomp(graph(U,'upper'), 'OutputForm', 'cell');
		
		if min(length(C{1}),length(C{2}))>1 | Msn(index(ind_i), index(ind_j)) == 1,

			% Level (function value) at which cluster is split
			ctr.Node{leaf}.split = mW;
			ctr.Node{leaf}.link = link;

			fmax = zeros(2,1);
			for j = 1:2,
				tmp = index(C{j});
				[fmax(j), mode] = max( fx(index(C{j})) );
				nd = struct('index', index(C{j}), 'level', mW, 'split', fmax(j), ...
					'mode', tmp(mode), 'depth',0, 'link', link);

				ctr = ctr.addnode(leaf, nd);
			end
			ctr.Node{leaf}.depth = min(fmax) -  mW;
			ctr.Node{leaf}.mode = [];
			break;

		% Isolated vertex with no local minima along connecting edge
		elseif max(length(C{1}),length(C{2})) > 1,
			id = 1;
			if length(C{2})==1,
				id = 2;
			end

			% remove isolated observation from adjacency matrix
			U= U(C{3-id}, C{3-id});
			index = index(C{3-id});

		% Final cluster has been identified
		elseif length(C{1})==1 & length(C{2})==1,
			if fx(index(C{1})) > fx(index(C{2})),
				ctr.Node{leaf}.mode = index(C{1});
			else
				ctr.Node{leaf}.mode = index(C{2});
			end
			l = l+1;
			break;

		else
			error('GPART: Error we should never be here\n');
		end
	end
end

%fprintf('partitioning finished\n');
%keyboard

%%%% Update Mode locations for internal tree nodes
ctr = update_modes(ctr,fx,1);

fcount = 0;
% PRUNE Cluster Tree
l = sort(ctr.findleaves(), 'descend');
checked = 0;
while length(l) > 1,

	p = ctr.getparent(l(1));
	s = ctr.getsiblings(l(1));
	pos = find(l==s);
	
	% mode of current leaf
	ml = ctr.Node{l(1)}.mode;
	% mode of sibling 
	ml = ctr.Node{l(1)}.mode;
	mr = ctr.Node{s}.mode;

	prune = false;
	% If sibling (S) is leaf or if mode of S has higher function value than mode of M
	if ~isempty(pos) | fx(mr) > fx(ml),
		ml = ctr.Node{l(1)}.mode;
		mr = ctr.Node{s}.mode;

		[prune, fc] = nomin(X(ml,:), X(mr,:), [fx(ml); fx(mr)], fhandle, pB); 
		fcount = fcount + fc;
	end

	if prune,
		% S is a leaf node
		if ~isempty(pos),
			% Update parent to become leaf
			ctr.Node{p}.depth = 0;

			ctr.Node{p}.split = max(ctr.Node{l(1)}.split, ctr.Node{s}.split);

			ctr.Node{p}.link = [0 0];

			ctr = ctr.chop(l(1));
			ctr = ctr.chop(s);

		% Sibling (s) is a subtree not a leaf
		else  
			% Remove leaf and its sibling effectively replacing split at parent
			% node with the split at the sibling of l(1)

			% Update parent: Parent node takes on properties of sibling
			ctr.Node{p}.split = ctr.Node{s}.split;
			ctr.Node{p}.depth = ctr.Node{s}.depth;
			ctr.Node{p}.link = ctr.Node{s}.link;

			% remove subtrees
			o = sort([l(1), s], 'descend');
			ctr = ctr.removenode(o(1));
			ctr = ctr.removenode(o(2));
			%keyboard
		end

		%%% DEBUGGING
		l = sort(ctr.findleaves(), 'descend');
		l = l(checked+1:end);
	else
		if ~isempty(pos),
			checked = checked+2;
			l = l([2:pos-1,pos+1:end]);
		else
			checked = checked+1;
			l = l([2:end]);
		end
	end
end

% Points in ROAs of objective function
l = ctr.findleaves();
idx = zeros(n,1);
for i=1:length(l),
	idx( ctr.Node{l(i)}.index ) = i;
end
noise = find(idx == 0);


% Noise Point Assignment 
outliers = sparse(n-length(noise),1);
l = [1];
pos = 1;
while length( setdiff(ctr.findleaves(), l) ) > 0,

	c = sort( ctr.getchildren(l(pos)) );
	if isempty(c),
		pos = pos+1;
	else
		% row vector containing indices of observations belonging to parent
		ip = ctr.Node{l(pos)}.index';
		% index of observations belonging to left child
		il = ctr.Node{c(1)}.index';
		% index of observations belonging to right child
		ir = ctr.Node{c(2)}.index';
		if length(ip) > length(il) + length(ir),
			% unallocated observations
			u = setdiff(ip, [ir, il]);
			out = [];
			for k = 1:length(u),
				% the edge with maximum weight is the one removed last
				w = max(nonzeros(A(u, [ir,il,out])));
				%if length(w) ~= 1, keyboard end

				% matrix is not symmetric but ties cannot be excluded 
				[i,j] = find(A(u, [ir,il,out]) == w);
				if isempty(i),
					fprintf('empty i\n');
					keyboard
				end
				%keyboard
				
				if j(1) <= length(ir),
					if Msn( u(i(1)), ir( j(1) ) ) == 1,
						outliers( u(i(1)) ) = u(i(1));
						out = [out, u(i(1))];
					else
						ir = [ir, u(i(1))];
					end

				elseif j(1) <= length(ir) + length(il),
					if Msn( u(i(1)), il( j(1) - length(ir) ) ) ==1,
						outliers( u(i(1)) ) = u(i(1));
						out = [out, u(i(1))];
					else
						il = [il, u(i(1))];
					end
				else
					%keyboard
					out = [out, u(i(1))];
					% assign to same cluster as previous outlier if minimum not found
					if Msn(u(i(1)), out( j(1) - length(ir) - length(il) )) ~= 1,
						outliers( u(i(1)) ) = outliers( out(j(1) - length(ir) - length(il)) );
					else
						outliers( u(i(1)) ) = u(i(1));
					end
				end
				u = [u(1:(i(1)-1)),u(i(1)+1:end)];
			end
			ctr.Node{c(1)}.index = sort(il)';
			ctr.Node{c(2)}.index = sort(ir)';
		end
		l = [l(1:pos-1), c, l(pos+1:end)];
	end
end

% Assignment of all observations to ROAs
l = ctr.findleaves();
K = length(l);
idx = zeros(n,1);
for i=1:K,
	idx(ctr.Node{l(i)}.index) = i;
end

%keyboard
uniq = unique(nonzeros(outliers));
for i = 1:length(uniq),
	idx(outliers == uniq(i)) = K+i;
end


if sum(idx==0) > 0,
	warning('Not all noise points are allocated to clusters!\n');
	%keyboard;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ctr = update_modes(ctr,fx,nid);

if isempty(ctr.Node{nid}.mode),
	children = ctr.getchildren(nid);
	if isempty(ctr.Node{children(1)}.mode),
		ctr = update_modes(ctr, fx, children(1));
	end
	if isempty(ctr.Node{children(2)}.mode),
		ctr = update_modes(ctr, fx, children(2));
	end

	if fx(ctr.Node{children(1)}.mode) > fx(ctr.Node{children(2)}.mode),
		ctr.Node{nid}.mode = ctr.Node{children(1)}.mode;
	else
		ctr.Node{nid}.mode = ctr.Node{children(2)}.mode;
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out, fcount] = nomin(x1,x2,fx,fhandle, np)

% find modes of corresponding clusters
fmin = min(fx);

t = linspace(0, 1, np + 2);
% the following ensures that points from in the midpoint and outwards are considered
order=ones(np,1);
order(1) = floor(np/2)+1;
order(2:2:end) = order(1) + cumsum(order(2:2:end));
order(3:2:end) = order(1) - cumsum(order(3:2:end));
order(np) = np+1;

fcount = 0;
out = true;
for l=1:np,
	f = fhandle( t(order(l))*x1 + (1. - t(order(l)) )*x2 );
	fcount = fcount + 1;
	if (f < fmin) & abs(f-fmin)>sqrt(eps),
		out = false;
		return;
	end
end
