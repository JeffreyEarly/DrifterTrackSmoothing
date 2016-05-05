function [m_x,Cm_x,B,Bq,tq,Wx] = bspline_fit_with_tension(t,x,sigma,t_knot,S,T,lambda,weight_function)
% drifter_fit_bspline    Find the maximum likelihood fit
%
% t         independent variable (time), length N
% x         observations at time t_i, length N
% sigma     error of observation x, scalar or length N vector
% t_knot    knot points for the splines
% S         degree of spline (e.g., 3 denotes cubic), scalar
% T         degree at which tension is applied (e.g., 2 denotes acc.)
% lambda    the tension parameter
% weight_function   used for iteratively reweighted least-squares
%
% The tension parameter lambda should be given as if it is (d-1)/(d*u^2).
% This function will then internally scale it by N/Q; d is the number of
% degrees of freedom and must be a value 1 or great (if it is 1, then
% you've set no tension). The value u corresponds to the rms-value of the
% T-th derivative.
% 
% B is a matrix containing the M B-splines, at the N locations t_i, for
% K=S+1 derivatives. The matrix is NxMxK.
%
% m_x       The coefficients for the model fit, length M.
% Cm_x      The error in the coefficients, size MxM
% B         The splines used for the observations
% Bq        The splines used for the quadrature (integration) grid.
% tq        The associated time points for the quadrature grid.
% Wx        The final weighting matrix used.

if (length(t) ~= length(x) )
   disp('The time series are not consistent lengths');
   return;
end

if (length(sigma) == 1)
   sigma = ones(size(x))*sigma; 
end

K = S+1;
B = bspline(t,t_knot,K);
X = squeeze(B(:,:,1));
N = size(X,1); % number of data points ;

% Now we need a quadrature (integration) grid that is finer
Q = 10*N; % number of points on the quadrature grid
tq = linspace(t(1),t(end),Q)';
Bq = bspline(tq,t_knot,K);

tension = zeros(S,1);
tension(T) = lambda;
gamma = tension*N/Q;

Wx = diag(1./(sigma.^2));

[m_x,Cm_x] = ComputeSolution( X, Bq, Wx, gamma, x );

error_x_previous = sigma;
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

