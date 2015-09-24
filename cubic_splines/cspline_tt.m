%---------------------------------------------------------
% cspline_tt
%---------------------------------------------------------
% calculate the 2nd derivative of cubic b-spline defined in Constable and
% Parker, 1997.
%
% The cubic spline is defined piecewise on the interval [-2,2] at knot
% points [-2,-1,0,1,2]. Everything outside of interval is zero.
%
% input: 	t	t value
% output:	x	spline 2nd derivative d2x/dt2 at x(t)
%
% Note that N. Teanby used t=(x-s)/k - (i - 2); where s is the starting
% point of the time series, k is the distance between knots and i is the
% index of the spline. This assumes that you want an extra spline on each
% end.
%
%---------------------------------------------------------
% N. Teanby	17/01/07
% N. Teanby	06/04/09	output curv changed from d2y/dt2 to d2y/dx2
% J. Early  18/09/2015  Modified to handle vector input.
% J. Early  18/09/2015  Removed the assumed normalization.
%---------------------------------------------------------
function [x] = cspline_tt(t)

x = zeros(size(t));
x = x + (-2 <= t & -1 > t).*6*(t+2);
x = x + (-1 <= t & 0 > t).*(-12 - 18*t);
x = x + (0 <= t & 1 > t).*(-12 + 18*t);
x = x + (1 <= t & 2 > t).*6*(2-t);
