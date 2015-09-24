%---------------------------------------------------------
% qspline_t - quintic spline first derivative
%---------------------------------------------------------
% Calculate the 1st derivative of the quintic b-spline defined in "A
% quintic B-spline finite-element method for solving the nonlinear
% Schrödinger equation" from B. Saka (2012).
%
% The quintic B-spline is defined piecewise on the interval [-3,3] with
% knot points [-3,-2,-1,0,1,2,3]. Everything outside that interval is zero.
%
% input: 	t	t value
% output:	x	spline 1st derivative dx/dt at x(t)
%
% Don't forget to renormalize the derivative by the knot spacing if it's
% not uniform in your problem.
%
%---------------------------------------------------------
% J. Early 14/09/2015
% J. Early 18/09/2015   Converted to vector format
%---------------------------------------------------------
function [x] = qspline_t(t)

x = zeros(size(t));
x = x + (-3 <= t & 3 > t).*( 5*(t+3).^4 );
x = x + (-2 <= t & 3 > t).*( -5*6*(t+2).^4 );
x = x + (-1 <= t & 3 > t).*( 5*15*(t+1).^4 );
x = x + ( 0 <= t & 3 > t).*( -5*20*(t+0).^4 );
x = x + ( 1 <= t & 3 > t).*(  5*15*(t-1).^4 );
x = x + ( 2 <= t & 3 > t).*( -5*6*(t-2).^4 );