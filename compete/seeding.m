function [idx, elapsed_time, seeds] = seeding(X,fhandle,radius,cecF)
% Reference:
% X. Li. Efficient Differential Evolution using Speciation for Multimodal Function Optimization
% GECCO '05 Proceedings of the 7th annual conference on Genetic and evolutionary computation, Pages 873-880, 2005

cec_rad = [0.01*ones(1,4), 0.5, 0.5, 0.2, 0.5, 0.2, 0.01*ones(1,11)];

if isempty(radius),
	radius = cec_rad(cecF);
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
ci = spalloc(n,1,n);
pos = 0;
[fs, order] = sort(fx,'descend');
for i = 1:n,
	
	%[fmax, j] = max(fx);
	% if minimum distance to existing seeds exceeds radius
	if i==1 || min( pdist2(X(order(i),:), X(ci(1:pos),:)) ) > radius,
		pos = pos + 1;
		ci(pos) = order(i);
	end
end
seeds = X(ci(1:pos),:);

idx = zeros(n,1);
for i=1:pos,
	rem = [1:(ci(i)-1), (ci(i)+1):n];
	index = find(pdist2( X(ci(i),:), X(rem,:) ) <= radius);
	index(index >= ci(i)) =  index(index >= ci(i)) + 1;
	idx(index) = i;
end
elapsed_time = toc(start_time);
