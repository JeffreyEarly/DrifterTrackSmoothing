%------------------------------------------------------------------------------
%
%     Constrained cbs curve fitting in tension
%     Nick Teanby 30/01/07
%
%------------------------------------------------------------------------------
%
%     Function to smooth a set of unevenly spaced x,y data
%	and output the cubic B spline parameters and covariance.
%
%     Allows constraints to be imposed on the smooth curve
%     by using the method of Lagrange multipliers.
%
%	Constraints [optional] can be on y, dy/dx, or d2y/dx2.
%
%	NB. curvature at ends of curve will be constrained to zero by default.
%
%     Tension is applied using a quadratic spring approximation as explained
%	in Teanby 2007.
%
%	  input
%       -----
%	x	float(n)		x data
%	dx	float(n)		x data errors
%
%	  input [optional]
%	  ----------------
%	x0	float(n0)		x constraints
%	y0	float(n0)		constraints [ y , dy/dx, d2y/dx2 ]
%	ctype	float(n0)		constraint type [0=y, 1=grad, 2=curvature]
%
%       output
%       ------
%	m	float(M)		spline parameters
%	cm	float(M,M)		covariance matrix of spline parameters
%
%------------------------------------------------------------------------------
function [m_x,Cm_x,B] = bspline_fit_no_tension(t,x,dx,S,t_knot,weight_function)

if (length(t) ~= length(x) )
    disp('The time series are not consistent lengths');
    return;
end

B = bspline(t,t_knot,S+1);
X = squeeze(B(:,:,1));
% N = size(X,1); % also size(X,1);
% M = size(X,2); % number of splines

Wx = diag(1./(dx.^2));

dbstop if warning
[m_x,Cm_x] = ComputeSolution( X, Wx, x );

error_x_previous = dx;
rel_error = 1.0;
repeats = 1;
while (rel_error > 0.01)
    dx2 = weight_function(X*m_x - x);
    
    Wx = diag(1./(dx2));
    
    [m_x,Cm_x] = ComputeSolution( X, Wx, x );
    
    rel_error = max( (dx2-error_x_previous)./dx2 );
    error_x_previous=dx2;
    repeats = repeats+1;
    
    if (repeats == 100)
        disp('Failed to converge after 100 iterations.');
        break;
    end
end

end


function [m_x, Cm_x] = ComputeSolution( X, Wx, x )
% X matrix:
% Rows are the N observations
% Columns are the M splines
% Wx is NxN
% x is Nx1

% set up inverse theory matrices
E_x = X'*Wx*X; % MxM

m_x = E_x\(X'*Wx*x);

Cm_x = inv(E_x);
end

