K = 4;
S = K-1;
t_data = linspace(0,10,11)';
t_data(4) = [];

N = length(t_data);
M = N-5;

DataDensityL = 1./diff(t_data);
DataDensityR = DataDensityL;
DataDensityL = [DataDensityL(1); DataDensityL];
DataDensityR = [DataDensityR; DataDensityR(end)];
DataDensity = min(DataDensityL,DataDensityR);
%DataDensity = (DataDensityL+DataDensityR)/2;

% DataDensity = ones(size(t_data));
% DataDensity(6) = 10;

t_knot = zeros(M,1);
for i=0:(M-1)
    tau_center = ((N-1)/(M-1))*i;
    half_width = ((N-1)/(M-1))*(K-2)/2;
    
    % shorten to the left boundary
    if (tau_center-half_width < 0)
        half_width = tau_center;
    end
    
    % shorten to the right boundary
    if (tau_center+half_width > (N-1))
        half_width = (N-1)-tau_center;
    end
    
    index_start = tau_center-half_width;
    index_end = tau_center+half_width;
    
    range = (floor(index_start):1:ceil(index_end))';
    weight = ones(size(range));
    
    if (mod(index_start,1) ~= 0 )
       weight(1) =  1-mod(index_start,1);
    end
        
    if (mod(index_end,1) ~= 0 )
       weight(end) =  mod(index_end,1);
    end
    
    weight = weight .* DataDensity(range+1);
    
    t_knot(i+1) = sum( t_data(range+1) .* weight ) / sum(weight);
end

figure
scatter(t_data,ones(size(t_data)))
xlim([-1 11])
ylim([0 2])
vlines(t_knot,'g--')

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