K = 6;
S = K-1;
t_data = linspace(0,10,11)';

N = length(t_data);
M = 11;

t_knot = zeros(M,1);
for i=1:M
    tau_center = (N/M)*i;
    half_width = (N/M)*(K-2)/2;
    if (tau_center-half_width < 1)
        half_width = tau_center-1;
    elseif (tau_center+half_width>N)
        half_width = N-tau_center;
    end
    
    index_start = tau_center-half_width;
    if (mod(index_start) ~= 0 )
       frac_start =  1-mod(index_end);
       index_start = floor(index_end);
    end
    
    index_end = tau_center+half_width;
    if (mod(index_end) ~= 0 )
       frac_end =  mod(index_end);
       index_end = ceil(index_end);
    end
    
    range = (tau_center-half_width):1:(tau_center+half_width);
    
    ceil(range(end))
    
    t_knot(i) = mean( t_data(range) );
end


return

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