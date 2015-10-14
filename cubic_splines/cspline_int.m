%---------------------------------------------------------
% cspline
%---------------------------------------------------------
% Calculate the value of cubic b-spline defined in Constable and Parker
% 1997
%
% The cubic spline is defined piecewise on the interval [-2,2] at knot
% points [-2,-1,0,1,2]. Everything outside of interval is zero.
%
% input: 	t	t value
% output:	x	spline value at x(t)
%
% Note that N. Teanby used t=(x-s)/k - (i - 2); where s is the starting
% point of the time series, k is the distance between knots and i is the
% index of the spline. This assumes that you want an extra spline on each
% end.
%
%---------------------------------------------------------
% N. Teanby	17/01/07
% J. Early  17/09/2015  Modified to handle vector input.
% J. Early  18/09/2015  Removed the assumed normalization.
%---------------------------------------------------------
function [x] = cspline_int(t)

x = zeros(size(t));
x = x + (-2 <= t & -1 > t).*(1/4)*(t+2).^4;
x = x + (-1 <= t & 0 > t).*(11/4 + 4*t - 2*t.^3 - (3/4)*t.^4);
x = x + (0 <= t & 1 > t).*(4*t - 2*t.^2 + (3/4)*t.^4);
x = x + (1 <= t & 2 > t).*(-(1/4)*(2-t).^4+(1/4));
