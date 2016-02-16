%---------------------------------------------------------
% qspline_tt - quintic spline second derivative
%---------------------------------------------------------
% Calculate the 2nd derivative of the quintic b-spline defined in "A
% quintic B-spline finite-element method for solving the nonlinear
% Schrödinger equation" from B. Saka (2012).
%
% The quintic B-spline is defined piecewise on the interval [-3,3] with
% knot points [-3,-2,-1,0,1,2,3]. Everything outside that interval is zero.
%
% input: 	t	t value
% output:	x	spline 2nd derivative dx2/dt2 at x(t)
%
% Don't forget to renormalize the derivative by the knot spacing if it's
% not uniform in your problem.
%
%---------------------------------------------------------
% J. Early 14/09/2015
% J. Early 18/09/2015   Converted to vector format
%---------------------------------------------------------
function [x] = qspline_tt(t)

x = zeros(size(t));
x = x + (-3 <= t & 3 > t).*( 4*5*(t+3).^3 );
x = x + (-2 <= t & 3 > t).*( -4*5*6*(t+2).^3 );
x = x + (-1 <= t & 3 > t).*( 4*5*15*(t+1).^3 );
x = x + ( 0 <= t & 3 > t).*( -4*5*20*(t+0).^3 );
x = x + ( 1 <= t & 3 > t).*(  4*5*15*(t-1).^3 );
x = x + ( 2 <= t & 3 > t).*( -4*5*6*(t-2).^3 );