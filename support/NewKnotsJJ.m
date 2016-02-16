% M is the number of new knot points that you want
% td is the location of the data points
% m_x and m_y are the coefficients
% tq is time on the quadrature (fine) grid
% Bq are the splines on the quadrature grid, at each derivative (K of them)
function [t_knot,t_knot_x,t_knot_y] = NewKnotsJJ( M, t_data, m_x, m_y, tq, Bq_x, Bq_y )

K = size(Bq_x,3);
N = length(t_data);

jx = squeeze(Bq_x(:,:,K))*m_x;
jy = squeeze(Bq_y(:,:,K))*m_y;

j = sqrt(jx.*jx+jy.*jy);

xi = cumsum(abs(j).^(1/K))*(tq(2)-tq(1));
j_data = interp1(tq,xi,t_data,'linear');
xi = cumsum(abs(jx).^(1/K))*(tq(2)-tq(1));
jx_data = interp1(tq,xi,t_data,'linear');
xi = cumsum(abs(jy).^(1/K))*(tq(2)-tq(1));
jy_data = interp1(tq,xi,t_data,'linear');


DataDensityL = 1./diff(t_data);
DataDensityR = DataDensityL;
DataDensityL = [DataDensityL(1); DataDensityL];
DataDensityR = [DataDensityR; DataDensityR(end)];
DataDensity = min(DataDensityL,DataDensityR);

t_knot = zeros(M,1);
t_knot(1) = t_data(1);
t_knot(end) = t_data(end);

t_knot_x = t_knot;
t_knot_y = t_knot;
for i=1:(M-2)
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
    
    range = range+1; % because Matlab indexes starting at 1
    t_knot(i+1) = sum( t_data(range) .* weight .* j_data(range) ) / sum(weight.*j_data(range));
    t_knot_x(i+1) = sum( t_data(range) .* weight .* jx_data(range) ) / sum(weight.*jx_data(range));
    t_knot_y(i+1) = sum( t_data(range) .* weight .* jy_data(range) ) / sum(weight.*jy_data(range));
end

return;











% Jq = squeeze(Bq_x(:,:,K));



t_knot = zeros(M,1);

for i=1:(M)
    a = i;
    b = i+K-1;
    if b>M
        b=M;
    elseif b>i+(i-1)
        b = i;
    end
    range = a:b;
    t_knot(i) = sum( t_data(range).*j_data(range) )/sum(j_data(range));
    t_knot_x(i) = sum( t_data(range).*jx_data(range) )/sum(jx_data(range));
    t_knot_y(i) = sum( t_data(range).*jy_data(range) )/sum(jy_data(range));
end
