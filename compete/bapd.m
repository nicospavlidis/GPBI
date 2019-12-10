function [idx, elapsed_time, dpeaks] = bapd(X,fhandle, eta)
% Reference: 
%
% Evolutionary Multiobjective Optimization-Based Multimodal Optimization:
% Fitness Landscape Approximation and Peak Detection.
%
% IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTATION, VOL. 22, NO. 5, OCTOBER 2018

if nargin < 3 || isempty(eta),
	eta = 0.1;
end

[n,d] = size(X);

global initial_flag;
initial_flag = 0;

fx = zeros(n,1);
for i=1:n,
	fx(i) = fhandle(X(i,:));
end
assert(sum(isnan(fx))==0, 'NaN in fx');
assert(sum(isinf(fx))==0, 'Infinite in fx');

start_time = tic;
ymin = min(fx);
ymax = max(fx);

ymin = eta*ymin + (1-eta)*ymax;
index = find(fx >= ymin);

peaks = {};
while length(index) > 2,
	%fprintf('length(index): %i\n', length(index));
	p = apd(X(index,:));

	% assign index in terms of original sequence
	for i=1:length(p),
		p{i} = index(p{i});
	end
	peaks = {peaks{:}, p{:}};

	ymin = 0.5*ymin + 0.5*ymax;
	while sum(fx >= ymin) == length(index) && (ymax-ymin)> 0.0001,
		ymin = 0.5*ymin + 0.5*ymax;
	end

	if ymax - ymin > 0.0001,
		index = find(fx >= ymin);
	else
		index = [];
	end
end

d = length(peaks);
iPeaks = zeros(1,d);
for i=1:length(peaks),
	[fmax, imax] = max( fx(peaks{i}) );
	% iPeaks(i): row index of maximiser of set peaks{i}
	iPeaks(i) = peaks{i}(imax);

	assert(fx(iPeaks(i)) - fmax == 0, 'Fmax not equal to corresponding fx entry');
end

%keyboard;
% indices of unique maximisers over all clusters (sorted in ascending order)
maxx = unique(iPeaks);
% X coordinates of maximisers (distinct peaks)
dpeaks = X(maxx, :);

% assign observations to clusters
idx = zeros(n,1);
for i = 1:length(maxx)
	% Find all clusters whose maximiser is the same as current maximiser
	id = find(iPeaks == maxx(i));
	% All other maximisers
	remain = [maxx(1:i-1), maxx(i+1:end)];

	% Find the largest cluster whose maximiser is the current maximiser
	% and which does not contain any other maximiser
	l1 = length(id);
	cluster = 0;
	sz = 0;
	% for each of the other subsets with the same maximiser
	for j = 1:length(id),
		% if the subset does not contain any other peak/ maximiser in
		% it and if it is the largest so far
		if isempty( intersect(remain, peaks{id(j)}) ) && length(peaks{id(j)}) > sz,
			sz = length( peaks{ id(j) } );
			cluster = id(j);
		end
	end

	if cluster ~= 0,
		idx( peaks{cluster} ) = i;
	else
		fprintf('BAPD: Strange behaviour\n');
	end
end
elapsed_time = toc(start_time);
%fprintf('Elapsed Time %1.5f\n', elapsed_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [peaks, sigmas] = apd(Xc)

n = size(Xc,1);
index = [1:n];

% l1 norm
D = pdist(Xc,'minkowski',1);
% find maximum nearest neighbour distance
%	sort columns in ascending order
[S,I] = sort(squareform(D), 1);

sigmas = zeros(ceil(n/2),1);
peaks = cell(ceil(n/2),1);
k=0;
while length(index) > 0,
	
	% maximum distance nearest neighbour
	% (since sigma can only decrease in each iteration it is impossible that 
	% the current "most distant NN" would be a member of any of the previous clusters
	% because if that was true then the point in index would also be a member of that cluster
	% since it is within a radius smaller than the previous sigma from another cluster member
	[sigma, i] = max( S(2,index) );

	% debugging
	for j=1:k,
		if sum(peaks{j} == index(i))~=0,
			fprintf('Starting point previously assigned\n');
			keyboard;
		end
		if sum(peaks{j} == I(2,index(i)))~=0,
			fprintf('NN assigned to previous cluster\n');
			keyboard;
		end
	end
	
	ci = spalloc(1,n,n);
	ci(1) = [index(i)];
	i = 1;
	% cluster length
	lc = 1;
	while i <= lc,
		% number of points (+1) whose distance to current observation <= s
		% It is impossible for points already allocated to a cluster to be encountered
		% again due to shrinking radius. Therefore consider all rows of column ci(i)
		np = sum( S(:, ci(i)) <= sigma );

		% Elements in first array not encountered in second
		pts = setdiff( I(1:np, ci(i)), ci(1:lc));
	
		if length(pts) > 0,
			% debugging
			for j=1:k,
				if ~isempty( intersect(peaks{j},pts) ),
					fprintf('This should never have occurred\n');
					keyboard;
				end
			end

			ci(lc+1:lc+length(pts)) = pts;
			lc = lc + length(pts);
		end
		i = i+1;
	end

	k = k+1;
	peaks{k} = ci(1:lc);
	sigmas(k) = sigma;
	index = setdiff(index, ci(1:lc));
end

peaks = {peaks{1:k}};
sigmas = sigmas(1:k);
%fprintf('Exiting APD\n');
