addpath('compete');
addpath('heatmaps');
addpath(genpath('cec13mmo'));

k=4;
b = 0.5*k; 
B = floor(b*1000);
C = 10;
fh = @(x)(niching_func(x,7));

T = table([],[],[],[],[],[],[],[],[],[],'VariableNames',{'method', 'iter','K', 'budget', 'NMI', 'aRand','purity',' numclusters','time','trueclusters'});
for iter=1:49,

	X = load( sprintf('DEpop/dl_fcn_7_strat_1_iter_%i.dat',iter) );
	[n,d] = size(X);
	% last column contains labels
	labels = X(:,d);
	X = X(:, 1:d-1);
	trueK = length(unique(labels));

	% NBC2
	[idx,elapsed_time, Wnbc2] = nbc2(X,fh);
	pn = performance(idx, labels);
	T = [T; {'nbc2',iter+1,0, 0, pn.NMI, pn.AdjRand, pn.Purity, max(idx), elapsed_time, trueK}];
	

	% Uncomment to produce heatmaps illustrating allocation to ROAs
	%cmp = {};
	%for i=1:max(idx), 
	%	cmp{i} = find(idx == i)'; 
	%end;
	%h = visual(7, X, Wnbc2, cmp, [],[0.25,0.25],[10.,0.82]);
	%fname = sprintf('heatmaps/VincentZoom_NBC2_iter_%i', iter);
	%print(strcat(fname,'.pdf'),'-dpdf','-bestfit');
	%set(gcf,'PaperPositionMode','auto'); print(fname,'-dpng','-r0')
	%close(h);

	% BAPD
	[idx,elapsed_time] = bapd(X,fh);
	pn = performance(idx, labels);
	T = [T; {'bapd',iter+1,0, 0, pn.NMI, pn.AdjRand, pn.Purity, max(idx), elapsed_time, trueK}];

	% SEEDS
	[idx,elapsed_time] = seeding(X,fh,[],7);
	pn = performance(idx, labels);
	T = [T; {'seed',iter+1,0, 0, pn.NMI, pn.AdjRand, pn.Purity, max(idx), elapsed_time, trueK}];

	[idx,noise,~,~,S,myclock] = gpbi(X, @(x)(niching_func(x,7)), k, B, C);

	p = performance(idx, labels);

	% Necessary for visualisation
	comp = {}; for i=1:max(idx), comp{i} = find(idx == i)'; end

	% Uncomment to visualise GPBI assignment of all points
	%h = visual(7, X, S, comp, [],[0.25,0.25],[10.,0.82]);
	%fname = sprintf('heatmaps/VincentZoom_all_k_%i_M_%i_iter_%i',k,B,iter);
	%print(strcat(fname,'.pdf'),'-dpdf','-bestfit');
	%set(gcf,'PaperPositionMode','auto'); print(fname,'-dpng','-r0')
	%close(h);

	if ~isempty(noise),
		labN = idx;
		labN(noise) = 0;
		l = unique(labN);
		estK = length(l)-1;

		% Necessary for visualisation
		cmp = {}; 
		for i=2:length(l),
			cmp{i-1} = find( labN == l(i) )';
		end
	else
		labN = idx;
		estK = max(idx);

		% Necessary for visualisation
		cmp = comp;
	end

	% Uncomment to visualise GPBI assignment of non-noise points
	%h = visual(7, X, S, cmp, [],[0.25,0.25],[10.,0.82]);
	%fname = sprintf('heatmaps/VincentZoom_noise_k_%i_M_%i_iter_%i',k,B,iter);
	%print(strcat(fname,'.pdf'),'-dpdf','-bestfit');
	%set(gcf,'PaperPositionMode','auto'); print(fname,'-dpng','-r0')
	%close(h);

	fprintf('Estimated Clusters including Outliers: %i\n', max(idx));
	fprintf('Estimated Clusters excluding Outliers: %i\n', estK);

	% graph partitioning for basin identification
	T = [T; {'gpbi',iter+1,k, B, p.NMI, p.AdjRand, p.Purity, estK, myclock, trueK}];
end
writetable(T,sprintf('VincentZoom_de_rand1_k%i_B%i.txt',k,b))
