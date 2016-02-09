function [m_x,Cm_x,B,Bq,tq] = bspline_fit_with_tension(t,x,dx,S,p, weight_function)
% drifter_fit_bspline    Find the maximum likelihood fit
%
% t         independent variable (time), length N
% x         observations at time t_i, length N
% dx        error of observation x, length N
% S         degree of spline (e.g., 3 denotes cubic), scalar
% t_knot    knot points for the splines
% p         the tension parameter
% weight_function   used for iteratively reweighted least-squares
%
% The tension parameter p should be given as if it is 1/u^2. This function
% will then internally scale it by N/Q;
% 
% B is a matrix containing the M B-splines, at the N locations t_i, for
% K=S+1 derivatives. The matrix is NxMxK.
%
% m_x       The coefficients for the model fit, length M.
% Cm_x      The error in the coefficients, size MxM
%
% 

if (length(t) ~= length(x) )
   disp('The time series are not consistent lengths');
   return;
end

% If we were only given one tension value, assume it's for the end.
if length(p) == 1
   pin = p;
   p = zeros(S,1);
   p(end) = pin;
end

if (length(p) < S)
    disp('The vector p must be of length S. You must sent a tension for 1st, 2nd, 3rd...S derivatives.');
    return;
end

K = S+1;

t_knot = NaturalKnotsForSpline( t, K );
B = bspline(t,t_knot,K);
X = squeeze(B(:,:,1));
N = size(X,1); % number of data points ;
M = size(X,2); % number of splines

% Now we need a quadrature (integration) grid that is finer
Q = 10*N; % number of points on the quadrature grid
tq = linspace(t(1),t(end),Q)';
Bq = bspline(tq,t_knot,K);

gamma = p*N/Q;


Wx = diag(1./(dx.^2));

[m_x,Cm_x] = ComputeSolution( X, Bq, Wx, gamma, x );

error_x_previous = dx;
rel_error = 1.0;
repeats = 1;
while (rel_error > 0.01)
    dx2 = weight_function(X*m_x - x);
    
    Wx = diag(1./(dx2));
    
    [m_x,Cm_x] = ComputeSolution( X, Bq, Wx, gamma, x );
    
    rel_error = max( (dx2-error_x_previous)./dx2 );
    error_x_previous=dx2;
    repeats = repeats+1;
    
    if (repeats == 100)
        disp('Failed to converge after 100 iterations.');
        break;
    end
end

end


function [m_x,Cm_x] = ComputeSolution( X, Bq, Wx, gamma, x )
% X matrix:
% Rows are the N observations
% Columns are the M splines
% Wx is NxN
% x is Nx1

E_x = X'*Wx*X; % MxM

for i=1:(size(Bq,3)-1)
    if (gamma(i) ~= 0.0)
        Xq = squeeze(Bq(:,:,i+1));
        T = gamma(i)*(Xq'*Xq);
        E_x = E_x + T;
    end
end

m_x = E_x\(X'*Wx*x);
Cm_x = inv(E_x);

end

