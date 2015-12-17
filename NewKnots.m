% M is the number of new knot points that you want
% m_x and m_y are the coefficients
% tq is time on the quadrature (fine) grid
% Bq are the splines on the quadrature grid, at each derivative
function t_knot = NewKnots( M, m_x, m_y, tq, Bq )

iDim = size(Bq,3);
Jq = squeeze(Bq(:,:,iDim));

jx2 = Jq*m_x;
jy2 = Jq*m_y;

j2 = sqrt(jx2.*jx2+jy2.*jy2);

xi = cumsum(abs(j2).^(1/(S+1)))*(tq(2)-tq(1));
xi_interp = linspace(xi(1),xi(end),M);
t_knot = interp1(xi,tq,xi_interp,'linear')';