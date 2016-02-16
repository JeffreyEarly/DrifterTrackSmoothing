%---------------------------------------------------------
% cspline_ttt
%---------------------------------------------------------
% calculate the 3rd derivative of cubic b-spline defined in Constable and
% Parker, 1997.
%
% The cubic spline is defined piecewise on the interval [-2,2] at knot
% points [-2,-1,0,1,2]. Everything outside of interval is zero.
%
% input: 	t	t value
% output:	x	spline 3rd derivative d3x/dt3 at x(t)
%
%---------------------------------------------------------
% J.Early	14/09/2015
% J. Early  18/09/2015  Modified to handle vector input.
%---------------------------------------------------------
function [x] = cspline_ttt(t)

x = zeros(size(t));
x = x + (-2 <= t & -1 > t)*6;
x = x + (-1 <= t & 0 > t)*(-18);
x = x + (0 <= t & 1 > t)*(+18);
x = x + (1 <= t & 2 > t)*(-6);