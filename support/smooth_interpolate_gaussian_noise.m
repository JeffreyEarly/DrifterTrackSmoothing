function [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t,x,y,sigma,S,T,lambda_x, lambda_y, DF)
% drifter_fit_bspline    Find the maximum likelihood fit
%
% t         independent variable (time), length N
% x,y       observations at time t_i, length N
% sigma     position error, std. dev. of Gaussian noise
% S         degree of spline (e.g., 3 denotes cubic splines, piecewise continuous accelerations), scalar
% T         degree at which tension is applied (e.g., 2 denotes acc.)
%
% The tension parameter lambda_x, lambda_y should be given as if it is (d-1)/(d*u^2). This function
% will then internally scale it by N/Q; d is the number of degrees of
% freedom and must be a value 1 or great (if it is 1, then you've set no
% tension). The value u corresponds to the rms-value of the T-th
% derivative.
% 
% B is a matrix containing the M B-splines, at the N locations t_i, for
% K=S+1 derivatives. The matrix is NxMxK.
%
% m_x       The coefficients for the model fit, length M.
% Cm_x      The error in the coefficients, size MxM
%
% 

if (length(t) ~= length(x) || length(t) ~= length(y) )
   disp('The time series are not consistent lengths.');
   return;
end

if (abs(diff(diff(t))) > 1e-6)
   disp('This function only works for evenly sampled data.');
   return;
end



if nargin < 9
    DF = 1;
end

K = S+1;
t_knot = NaturalKnotsForSpline( t, K, DF );
B = bspline(t,t_knot,K);
X = squeeze(B(:,:,1));
N = size(X,1); % number of data points ;

% Now we need a quadrature (integration) grid that is finer
Q = 10*N; % number of points on the quadrature grid
tq = linspace(t(1),t(end),Q)';
Bq = bspline(tq,t_knot,K);

tension_x = zeros(S,1);
tension_y = zeros(S,1);

tension_x(T) = lambda_x;
tension_y(T) = lambda_y;
gamma_x = tension_x*N/Q;
gamma_y = tension_y*N/Q;


[m_x,m_y,Cm_x,Cm_y] = ComputeSolution( X, Bq, sigma, gamma_x, gamma_y, x, y );


end

function totalError = TotalPositionErrorRatio(  X, Bq, sigma, x, y, a, T, N, Q )

a = 10^(a);
tension = zeros((size(Bq,3)-1),1);
tension(T) = 1/a^2;
gamma = tension*N/Q;
[m_x,m_y,~,~] = ComputeSolution( X, Bq, sigma, gamma, x, y );

x_out = X*m_x - x; y_out = X*m_y - y;

lambda = mean( [x_out.*x_out; y_out.*y_out] )/sigma^2;

shouldDisplay = 1;
if shouldDisplay == 1
    fprintf('\t(acceleration, lambda) = (%g, %f)\n', a, lambda);
end

totalError = abs(lambda-1.0);
%     totalError = abs(log10(lambda));

end

function totalError = TotalPowerRatio(  X, Bq, sigma, x, y, a, T, D, N, Q, observed_power, noise_power )

a = 10^(a);
tension = zeros((size(Bq,3)-1),1);
tension(T) = 1/a^2;
gamma = tension*N/Q;
[m_x,m_y,~,~] = ComputeSolution( X, Bq, sigma, gamma, x, y );

x_out = X*m_x; y_out = X*m_y;

tensioned_power = mean( (D*x_out).^2 + (D*y_out).^2 );

lambda = (observed_power - tensioned_power)/noise_power;

shouldDisplay = 1;
if shouldDisplay == 1
    fprintf('\t(acceleration, lambda) = (%g, %f)\n', a, lambda);
end

totalError = abs(lambda-1.0);
%     totalError = abs(log10(lambda));

end

function [m_x,m_y,Cm_x,Cm_y] = ComputeSolution( X, Bq, sigma, gamma_x, gamma_y, x, y )
% X matrix:
% Rows are the N observations
% Columns are the M splines
% Wx is NxN
% x is Nx1

E_x = X'*X/(sigma*sigma); % MxM
E_y = X'*X/(sigma*sigma); % MxM

for i=1:(size(Bq,3)-1)
    if (gamma_x(i) ~= 0.0)
        Xq = squeeze(Bq(:,:,i+1));
        T = gamma_x(i)*(Xq'*Xq);
        E_x = E_x + T;
    end
    if (gamma_y(i) ~= 0.0)
        Xq = squeeze(Bq(:,:,i+1));
        T = gamma_y(i)*(Xq'*Xq);
        E_y = E_y + T;
    end
end


m_x = E_x\(X'*x/(sigma*sigma));
m_y = E_y\(X'*y/(sigma*sigma));

Cm_x = inv(E_x);
Cm_y = inv(E_y);
end

