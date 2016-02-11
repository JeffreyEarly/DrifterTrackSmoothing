function [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = bspline_bivariate_fit_with_tension(t,x,y,dx,dy,S,p, weight_function, t_knot)
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

if (length(t) ~= length(x) || length(t) ~= length(y) )
   disp('The time series are not consistent lengths');
   return;
end

% If we were only given one tension value, assume it's for the end.
if length(p) == 1
   pin = p;
   p = zeros(S-1,1);
   p(end) = pin;
end

if (length(p) < S-1)
    disp('The vector p must be of length S-1. You must sent a tension for 1st, 2nd, 3rd...S-1 derivatives.');
    return;
end

K = S+1;

if nargin < 9
    t_knot = NaturalKnotsForSpline( t, K );
end
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
Wy = diag(1./(dy.^2));

dbstop if warning
[m_x,m_y,Cm_x,Cm_y] = ComputeSolution( X, Bq, Wx, Wy, gamma, x, y );

error_x_previous = dx;
error_y_previous = dy;
rel_error = 1.0;
repeats = 1;
while (rel_error > 0.01)
    dx2 = weight_function(X*m_x - x);
    dy2 = weight_function(X*m_y - y);
    
    Wx = diag(1./(dx2));
    Wy = diag(1./(dy2));
    
    [m_x,m_y,Cm_x,Cm_y] = ComputeSolution( X, Bq, Wx, Wy, gamma, x, y );
    
    rel_error = max(max( (dx2-error_x_previous)./dx2 ),max( (dy2-error_y_previous)./dy2 ));
    error_x_previous=dx2;
    error_y_previous=dy2;
    repeats = repeats+1;
    
    if (repeats == 100)
        disp('Failed to converge after 100 iterations.');
        break;
    end
end

dx=X*m_x - x;
dy=X*m_y - y;
a=dx'*Wx*dx/N;
b=dy'*Wy*dy/N;
chi = 2/(a+b);
sigma2 = sum(dx.*dx + dy.*dy)/(2*N);
w_rms = median(sqrt(1./[diag(Wx);diag(Wy)]));
fprintf('sum(dx/sigma)^2=%f, sum(dy/sigma)^2=%f, chi=%f, sigma = %f, w_rms=%f\n',a,b,chi,sqrt(sigma2),w_rms);


end


function [m_x,m_y,Cm_x,Cm_y] = ComputeSolution( X, Bq, Wx, Wy, gamma, x, y )
% X matrix:
% Rows are the N observations
% Columns are the M splines
% Wx is NxN
% x is Nx1

E_x = X'*Wx*X; % MxM
E_y = X'*Wy*X; % MxM

for i=1:(size(Bq,3)-1)
    if (gamma(i) ~= 0.0)
        Xq = squeeze(Bq(:,:,i+1));
        T = gamma(i)*(Xq'*Xq);
        E_x = E_x + T;
        E_y = E_y + T;
    end
end


m_x = E_x\(X'*Wx*x);
m_y = E_y\(X'*Wy*y);

Cm_x = inv(E_x);
Cm_y = inv(E_y);
end

