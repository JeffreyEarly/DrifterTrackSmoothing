% M is the number of new knot points that you want
% m_x and m_y are the coefficients
% tq is time on the quadrature (fine) grid
% Bq are the splines on the quadrature grid, at each derivative
function [t_knot,t_knot_x,t_knot_y] = NewKnots( M, m_x, m_y, tq, Bq_x, Bq_y )

K = size(Bq_x,3);
% Jq = squeeze(Bq_x(:,:,K));

jx = squeeze(Bq_x(:,:,K))*m_x;
jy = squeeze(Bq_y(:,:,K))*m_y;

j = sqrt(jx.*jx+jy.*jy);

xi = cumsum(abs(j).^(1/K))*(tq(2)-tq(1));
xi_interp = linspace(xi(1),xi(end),M);
t_knot = interp1(xi,tq,xi_interp,'linear')';

xi = cumsum(abs(jx).^(1/K))*(tq(2)-tq(1));
xi_interp = linspace(xi(1),xi(end),M);
t_knot_x = interp1(xi,tq,xi_interp,'linear')';

xi = cumsum(abs(jy).^(1/K))*(tq(2)-tq(1));
xi_interp = linspace(xi(1),xi(end),M);
t_knot_y = interp1(xi,tq,xi_interp,'linear')';