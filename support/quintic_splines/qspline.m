%---------------------------------------------------------
% qspline - quintic spline
%---------------------------------------------------------
% Calculate the quintic b-spline defined in "A quintic B-spline
% finite-element method for solving the nonlinear Schrödinger equation"
% from B. Saka (2012).
%
% The quintic B-spline is defined piecewise on the interval [-3,3] with
% knot points [-3,-2,-1,0,1,2,3]. Everything outside that interval is zero.
%
% input: 	t	t value
% output:	x	spline at x(t)
%
%---------------------------------------------------------
% J. Early 14/09/2015
% J. Early 18/09/2015   Converted to vector format
%---------------------------------------------------------
function [x] = qspline(t)

x = zeros(size(t));
x = x + (-3 <= t & 3 > t).*( (t+3).^5 );
x = x + (-2 <= t & 3 > t).*( -6*(t+2).^5 );
x = x + (-1 <= t & 3 > t).*( 15*(t+1).^5 );
x = x + ( 0 <= t & 3 > t).*( -20*(t+0).^5 );
x = x + ( 1 <= t & 3 > t).*(  15*(t-1).^5 );
x = x + ( 2 <= t & 3 > t).*( -6*(t-2).^5 );
