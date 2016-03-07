function [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t,x,y,sigma,S,T,a_in)
% drifter_fit_bspline    Find the maximum likelihood fit
%
% t         independent variable (time), length N
% x,y       observations at time t_i, length N
% sigma     position error, std. dev. of Gaussian noise
% S         degree of spline (e.g., 3 denotes cubic splines, piecewise continuous accelerations), scalar
% T         degree at which tension is applied (e.g., 2 denotes acc.)
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
   disp('The time series are not consistent lengths.');
   return;
end

if (abs(diff(diff(t))) > 1e-6)
   disp('This function only works for evenly sampled data.');
   return;
end

dt = t(2)-t(1);
D = FiniteDifferenceMatrixNoBoundary(T,t,1);
u_signal = D*x;
v_signal = D*y;
observed_rms_power = sqrt( mean( u_signal.*u_signal + v_signal.*v_signal));
fprintf('The rms value of the %d-th derivative of the signal is %g\n',T, observed_rms_power)

% The extra sqrt 2 comes from combining both x and y directions
noise_coefficients = [2;6;20;70;252;924;3432];
noise_rms_power = sqrt(2) * sqrt(noise_coefficients(T))*sigma/dt^T;
fprintf('The rms value of the of the noise is %g\n', noise_rms_power)

if (noise_rms_power > observed_rms_power)
   disp('The total assumed noise variance is greater than the observed signal.');
   a = observed_rms_power;
else
   a = sqrt( observed_rms_power.*observed_rms_power - noise_rms_power.*noise_rms_power )/sqrt(2);
   fprintf('The deduced value of the tension is %g\n', a);
end

tension = zeros(S,1);
tension(T) = 1/a^2;

% Don't want to base this on the highest derivative, honestly
% Gamma = noise_rms_power/a;
% if round(Gamma/2) > 1
%     DF = round(Gamma/2);
%     fprintf('Reducing the total number of knot points. Setting DF=%d.\n',DF);
% else
%     DF = 1;
% end

DF = 1;

K = S+1;
t_knot = NaturalKnotsForSpline( t, K, DF );
B = bspline(t,t_knot,K);
X = squeeze(B(:,:,1));
N = size(X,1); % number of data points ;

% Now we need a quadrature (integration) grid that is finer
Q = 10*N; % number of points on the quadrature grid
tq = linspace(t(1),t(end),Q)';
Bq = bspline(tq,t_knot,K);

dbstop if warning

if nargin < 7
    
    fprintf('Seaching for a maximum tension parameter to preserve total variance...\n')
    errorFunction = @(a) TotalPositionErrorRatio(  X, Bq, sigma, x, y, a, T, N, Q );
    optimalAcceleration = fminsearch( errorFunction, log10(a), optimset('TolX', 0.01, 'TolFun', 0.01) );
    a = 10^(optimalAcceleration(1));
    fprintf('Optimal acceleration tension is %g\n', a );  
    
    % fprintf('Seaching for a new optimal tension parameter to preserve total variance...\n')
    % errorFunction = @(a) TotalPowerRatio(  X, Bq, sigma, x, y, a, T, D, N, Q, observed_rms_power^2, noise_rms_power^2 );
    % optimalAcceleration = fminsearch( errorFunction, log10(a), optimset('TolX', 0.01, 'TolFun', 0.01) );
    % fprintf('Optimal acceleration tension is %g\n', 10^(optimalAcceleration(1)) );
    % a = 10^(optimalAcceleration(1));
    
    %a = 0.2;
    % a = 0.0021;
    % a = 0.5*4.175e-04;
else
    a = a_in;
end

tension(T) = 1/a^2;
gamma = tension*N/Q;


[m_x,m_y,Cm_x,Cm_y] = ComputeSolution( X, Bq, sigma, gamma, x, y );


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

function [m_x,m_y,Cm_x,Cm_y] = ComputeSolution( X, Bq, sigma, gamma, x, y )
% X matrix:
% Rows are the N observations
% Columns are the M splines
% Wx is NxN
% x is Nx1

E_x = X'*X/(sigma*sigma); % MxM
E_y = X'*X/(sigma*sigma); % MxM

for i=1:(size(Bq,3)-1)
    if (gamma(i) ~= 0.0)
        Xq = squeeze(Bq(:,:,i+1));
        T = gamma(i)*(Xq'*Xq);
        E_x = E_x + T;
        E_y = E_y + T;
    end
end


m_x = E_x\(X'*x/(sigma*sigma));
m_y = E_y\(X'*y/(sigma*sigma));

Cm_x = inv(E_x);
Cm_y = inv(E_y);
end

