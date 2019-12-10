function labels = labelVincent(X)
%
% Function used to label any point of the 2-dimensional Vincent function in the
% range illustrated in Figure 2 of the paper

n = size(X,1);
labels = zeros(n,1);
for i = 1:n,
	if X(i,2) < 0.46,
		if X(i,1) < 0.45,
			labels(i) = 1;
		elseif X(i,1) < 0.85,
			labels(i) = 2;
		elseif X(i,1) < 1.60,
			labels(i) = 3;
		elseif X(i,1) < 3.0,
			labels(i) = 4;
		elseif X(i,1) < 5.75,
			labels(i) = 5;
		else
			labels(i) = 6;
		end
	else
		if X(i,1) < 0.45,
			labels(i) = 7;
		elseif X(i,1) < 0.85,
			labels(i) = 8;
		elseif X(i,1) < 1.60,
			labels(i) = 9;
		elseif X(i,1) < 3.0,
			labels(i) = 10;
		elseif X(i,1) < 5.75,
			labels(i) = 11;
		else
			labels(i) = 12;
		end
	end
end
