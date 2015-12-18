K = 4;
S = K-1;
t_knot = linspace(0,10,11)';

t_data = t_knot;

% rows are t_data, columns are number of splines
B = bspline(t_data,t_knot,K);
X = squeeze(B(:,:,1));

t_knot = [repmat(t_knot(1),S,1); t_knot; repmat(t_knot(end),S,1)];
t_knot2 = t_knot;
for i=1:(length(t_data)-K+2)
    range = (i+1-1):(i+K-1-1);
    t_knot2(K+i) = mean(t_data(range));
end

% each data point has K-1 nonzero splines, the values of which add to one.
% t_knot'*X