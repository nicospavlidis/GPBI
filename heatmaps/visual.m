function h = visual(cecFun,X,W,comp, optima,xlow,xup, Msn);
%
% cecFun: number of cec 2013 test function (see technical report)
% X: Matrix containing locations of minimisers
% W: Adjacency matrix
% comp: connected components cell array

if nargin<4,
	comp = [];
end

%addpath('cec2013');
global initial_flag;
initial_flag=0;
fhandle = @(x)( niching_func(x,cecFun) ); 


if nargin < 7 | isempty(xup),
	xup = get_ub(cecFun);
end
if nargin < 6 | isempty(xlow),
	xlow = get_lb(cecFun);
end

if nargin < 8 | isempty(Msn),
	Msn = sparse(size(W,1), size(W,2));
end

x = linspace(xlow(1), xup(1), 500);
y = linspace(xlow(2), xup(2), 500);
[X1,Y1] = meshgrid(x,y);

% the files produced 
	
if nargin < 6 & isfile( sprintf('heatmaps/%i_Z.dat', cecFun) ),
	%fprintf('Data load: full range\n');
	Z1 = load(sprintf('heatmaps/%i_Z.dat', cecFun));
	
elseif nargin >=6 & isfile( sprintf('heatmaps/%i_Z_zoom_%1.2f_%1.2f_%1.2f_%1.2f.dat', cecFun, ...
	xlow(1),xup(1),xlow(2),xup(2)) ),
	%fprintf('Data load: custom range\n');
	Z1 = load(sprintf('heatmaps/%i_Z_zoom_%1.2f_%1.2f_%1.2f_%1.2f.dat', ...
			cecFun, xlow(1),xup(1),xlow(2),xup(2)));
else
	Z1 = zeros(length(y), length(x));
	for i=1:length(y),
		for j=1:length(x),
			%Z1(i,j) = fhandle( [X1(i,j); Y1(i,j)] );
			Z1(i,j) = niching_func( [X1(i,j) Y1(i,j)], cecFun);
		end
	end
	if nargin<6,
		dlmwrite(sprintf('heatmaps/%i_Z.dat',cecFun), Z1,'delimiter',',');
	else
		dlmwrite(sprintf('heatmaps/%i_Z_zoom_%1.2f_%1.2f_%1.2f_%1.2f.dat', ... 
			cecFun, xlow(1),xup(1),xlow(2),xup(2)), Z1,'delimiter',',');
	end
end


%%%% VISUALISATION
% Plot minimisers and graph induced by W
h = figure('rend','painters','pos',[0 0 1600 400]);
hold on;
%%surfc(X1,Y1,Z1, 'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
imagesc(x,y,Z1);
colormap(jet(80));
%colorbar

% if a cell array containing connected components is provided as input
if ~isempty(comp),
	assigned = [];
	col = lines(length(comp));
	%col = prism(length(comp));
	for i=1:length(comp);
		Y = X(comp{i},:);
		if length(comp{i})>1,
			cl = col(i,:);
			%cl = 'w';
		else
			%cl = 'w';
			cl = 'k';
		end
		%plot(Y(:,1), Y(:,2), 'o', 'MarkerSize',3,'Color',col(1,:),'MarkerFaceColor', col(1,:));

		Wcomp = W(comp{i},comp{i});
		Mt = Msn(comp{i},comp{i});
		for j=1:size(Y,1),
			index = find(Wcomp(:,j) ~= 0);
			for k=1:length(index),
				Y1 = [Y(j,:); Y(index(k),:)];
				if Mt(j,index(k)) == 2, % no fcn evaluation
					plot(Y1(:,1), Y1(:,2), '-', 'Color', 'k');
				elseif Mt(j,index(k)) == 3, % no local minimum found
					plot(Y1(:,1), Y1(:,2), '-', 'Color', 'b');
				elseif Mt(j,index(k)) == 1, % local minimum 
					plot(Y1(:,1), Y1(:,2), ':', 'Color', 'r');
				else
					% M not specified
					plot(Y1(:,1), Y1(:,2), '-', 'Color', 'k','LineWidth',0.5);
					%plot(Y1(:,1), Y1(:,2), ':', 'Color', col(1,:),'LineWidth',0.4);
				end
			end
		end	   
		plot(Y(:,1), Y(:,2), 'o', 'MarkerSize',4,'Color', 'k','MarkerFaceColor',cl);
		assigned = [assigned, comp{i}];
	end

	% plot unassigned (noise) points
	plot(X(setdiff([1:size(X,1)]',assigned),1), X(setdiff([1:size(X,1)]',assigned),2), 'o', ...
		'MarkerSize',3,'Color','k','MarkerFaceColor','w');
	
else
	plot(X(:,1), X(:,2), 'o', 'MarkerSize',3)
	for i=1:size(X,1),
		index = find(W(:,i) ~= 0);
		for j=1:length(index),

			Y = [X(i,:); X(index(j),:)];
			if Msn(i,index(j)) == 2, % no fcn evaluation
				plot(Y(:,1), Y(:,2), ':', 'Color', 'k');
			elseif Msn(i,index(j)) == 3, % no local minimum found
				plot(Y(:,1), Y(:,2), ':', 'Color', 'c');
			elseif Msn(i,index(j)) == 1, % local minimum 
				plot(Y(:,1), Y(:,2), ':', 'Color', 'r');
			else
				% M not specified
				plot(Y(:,1), Y(:,2), '-', 'Color', 'b');
			end
			%plot(Y(:,1), Y(:,2), 'b:');
		end
	end

end

if nargin> 4 & ~isempty(optima),
	plot(optima(:,1), optima(:,2), '+', 'MarkerSize',3, 'Color','k');
end

hold off;
axis([xlow(1), xup(1), xlow(2), xup(2)])
drawnow
