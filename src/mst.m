function tr = mst(A)

n = size(A,1);
tr = spalloc(n,n,2*(n-1));

conn_nodes = 1;        % nodes part of the min-span-tree
rem_nodes = [2:n];     % remaining nodes

for k=1:n-1,
	
	link = min( nonzeros(A(conn_nodes, rem_nodes)) );
	if isempty(link),
		error('Adjacency matrix is disconnected');
	end
	[ind_i, ind_j] = find(A(conn_nodes, rem_nodes) == link);
	%[ind_i,ind_j] = ind2sub([k,n-k],ind(1));

	i = conn_nodes(ind_i(1));
	j = rem_nodes(ind_j(1));

	tr(i,j) = link;
	tr(j,i) = link;
	
	conn_nodes = [conn_nodes, j];
	rem_nodes = rem_nodes([1:(ind_j(1)-1), (ind_j(1)+1):end]);
end
