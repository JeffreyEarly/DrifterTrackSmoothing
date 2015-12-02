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
function [x] = cspline(t,varargin)

if isempty(varargin)
    k = [-2,-1,0,1,2];
elseif length(varargin) == 1
    k = varargin{1};
    if length(k) ~= 5
        disp('You must provide exactly 5 knot points for a cubic spline');
        return;
    end
end

x = zeros(size(t));
x = x + (k(1) <= t & k(2) > t).*(t+2).^3;
x = x + (k(2) <= t & k(3) > t).*(4 - 6*t.^2 - 3*t.^3);
x = x + (k(3) <= t & k(4) > t).*(4 - 6*t.^2 + 3*t.^3);
x = x + (k(4) <= t & k(5) > t).*((2-t).^3);
