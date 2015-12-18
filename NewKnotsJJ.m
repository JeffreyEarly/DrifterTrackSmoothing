% M is the number of new knot points that you want
% td is the location of the data points
% m_x and m_y are the coefficients
% tq is time on the quadrature (fine) grid
% Bq are the splines on the quadrature grid, at each derivative (K of them)
function [t_knot,t_knot_x,t_knot_y] = NewKnotsJJ( M, t_data, m_x, m_y, tq, Bq_x, Bq_y )

K = size(Bq_x,3);
% Jq = squeeze(Bq_x(:,:,K));

jx = squeeze(Bq_x(:,:,K))*m_x;
jy = squeeze(Bq_y(:,:,K))*m_y;

j = sqrt(jx.*jx+jy.*jy);

xi = cumsum(abs(j).^(1/K))*(tq(2)-tq(1));
j_data = interp1(tq,xi,t_data,'linear');
xi = cumsum(abs(jx).^(1/K))*(tq(2)-tq(1));
jx_data = interp1(tq,xi,t_data,'linear');
xi = cumsum(abs(jy).^(1/K))*(tq(2)-tq(1));
jy_data = interp1(tq,xi,t_data,'linear');

t_knot = zeros(M,1);
t_knot_x = zeros(M,1);
t_knot_y = zeros(M,1);
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
