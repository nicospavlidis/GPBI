function W = cgraph(X,y,k)
%Connectivity Graph
%function W = cgraph(X,y,k)
%
%Inputs:
% X: Matrix whose rows correspond to points in domain of f
% y: Function values
% k: Number of nearest neighbours to which each point is connected

[n,d] = size(X);

W = sparse(n,n);
for i=1:n,
	index = [1:i-1,i+1:n]';
	% identify kNNs
	d = pdist2(X(i,:), X(index,:));

	% dk := distance to k-NNs
	% nn := id of k-NNs
	[dk, nn] = mink(d,k);

	% Correct for location on k-NNs wrt location of i
	nn(nn >= i) = nn(nn >= i) + 1;
	
	% if none of kNNs are better
	if sum(y(nn) > y(i)) == 0,
		% index of better solutions
		bt = find(y(index) > y(i));
		% nearest-better	
		[db, nb] = min( d(bt) );
		nb = index(bt(nb));

		%if (y(nb) <= y(i)) | abs(norm(X(i,:) - X(nb,:),2) - db) > sqrt(eps),
		%	fprintf('Failure to get Nearest Better correctly\n');
		%	keyboard
		%end
		nn = [nn, nb];
	end
	
	W(nn,i) = 1;
end
% Make affinity matrix symmetric
W = max(W,W');
