%---------------------------------------------------------
% cspline_t
%---------------------------------------------------------
% calculate the first derivative of cubic b-spline defined in Constable and
% Parker, 1997.
%
% The cubic spline is defined piecewise on the interval [-2,2] at knot
% points [-2,-1,0,1,2]. Everything outside of interval is zero.
%
% input: 	t	t value
% output:	x	spline 1st derivative  dx/dt at x(t)
%
% Note that N. Teanby used t=(x-s)/k - (i - 2); where s is the starting
% point of the time series, k is the distance between knots and i is the
% index of the spline. This assumes that you want an extra spline on each
% end.
%
% Also, don't forget that the derivative must be normalized by the knot
% spacing, e.g. x = x/k;
%---------------------------------------------------------
% N. Teanby	17/01/07
% N. Teanby	06/04/09	output grad changed from dy/dt to dy/dx
% J. Early  17/09/2015  Modified to handle vector input.
% J. Early  18/09/2015  Removed the assumed normalization.
%---------------------------------------------------------
function [x] = cspline_t(t)

x = zeros(size(t));
x = x + (-2 <= t & -1 > t).*(3*(t+2).^2);
x = x + (-1 <= t & 0 > t).*(-12*t - 9*t.^2);
x = x + (0 <= t & 1 > t).*(-12*t + 9*t.^2);
x = x + (1 <= t & 2 > t).*(-3*(2-t).^2);
