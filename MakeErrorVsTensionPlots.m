addpath('./support');
% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

load('sample_data/SyntheticTrajectories.mat')


D = FiniteDifferenceMatrixNoBoundary(1,t,1);
D2 = FiniteDifferenceMatrixNoBoundary(2,t,1);
D3 = FiniteDifferenceMatrixNoBoundary(3,t,1);
dt = t(2)-t(1);
cv = D*(x + sqrt(-1)*y);
ca = D2*(x + sqrt(-1)*y);
cj = D3*(x + sqrt(-1)*y);

u_true = std(cv)/sqrt(2);
a_true = std(ca)/sqrt(2);
j_true = std(cj)/sqrt(2);

fprintf('Measured std-u=%g, measured std-a=%g, measured std-j=%g\n', u_true, a_true, j_true);

stride = 10;
indices = 1:stride:floor(length(t)/1);
indicesAll = 1:max(indices);
x_obs = x(indices) + epsilon_x(indices);
y_obs = y(indices) + epsilon_y(indices);
t_obs = t(indices);
sigma = position_error;
S = 3;
T = 1;
K = S+1;

% t_knot = NaturalKnotsForSpline( t, K, 1 );
% B_all = bspline(t,t_knot,K);

DF = 1;

u_vals = 10.^(linspace(-1.5,0,15))';
rms_error_T1 = zeros(size(u_vals));
chi2_est_T1 = zeros(size(u_vals));
final_variance_T1 = zeros(size(u_vals));
for i = 1:length(u_vals)
   [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t_obs,x_obs,y_obs,sigma,S,T,u_vals(i));
   X = squeeze(B(:,:,1));
   t_knot = NaturalKnotsForSpline( t_obs, K, DF );
   B = bspline(t(indicesAll),t_knot,K);
   X1 = squeeze(B(:,:,1));
   U1 = squeeze(B(:,:,T+1));
   rms_error_T1(i) = std( (x(indicesAll)-X1*m_x) + sqrt(-1)*(y(indicesAll)-X1*m_y) );
   chi2_est_T1(i) = std( (x_obs-X*m_x) + sqrt(-1)*(y_obs-X*m_y) );
   final_variance_T1(i) = std( U1*m_x + sqrt(-1)*U1*m_y )/sqrt(2);
end

T = 2;
a_vals = 10.^(linspace(-4,-2,15))';
rms_error_T2 = zeros(size(a_vals));
chi2_est_T2 = zeros(size(a_vals));
final_variance_T2 = zeros(size(a_vals));
for i = 1:length(a_vals)
   [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t_obs,x_obs,y_obs,sigma,S,T,a_vals(i));
   X = squeeze(B(:,:,1));
   t_knot = NaturalKnotsForSpline( t_obs, K, DF );
   B = bspline(t(indicesAll),t_knot,K);
   X1 = squeeze(B(:,:,1));
   U1 = squeeze(B(:,:,T+1));
   rms_error_T2(i) = std( (x(indicesAll)-X1*m_x) + sqrt(-1)*(y(indicesAll)-X1*m_y) );
   chi2_est_T2(i) = std( (x_obs-X*m_x) + sqrt(-1)*(y_obs-X*m_y) );
   final_variance_T2(i) = std( U1*m_x + sqrt(-1)*U1*m_y )/sqrt(2);
end

T = 3;
j_vals = 10.^(linspace(-7,-4,15))';
rms_error_T3 = zeros(size(j_vals));
chi2_est_T3 = zeros(size(j_vals));
final_variance_T3 = zeros(size(j_vals));
for i = 1:length(j_vals)
   [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t_obs,x_obs,y_obs,sigma,S,T,j_vals(i));
   X = squeeze(B(:,:,1));
   t_knot = NaturalKnotsForSpline( t_obs, K, DF );
   B = bspline(t(indicesAll),t_knot,K);
   X1 = squeeze(B(:,:,1));
   U1 = squeeze(B(:,:,T+1));
   rms_error_T3(i) = std( (x(indicesAll)-X1*m_x) + sqrt(-1)*(y(indicesAll)-X1*m_y) );
   chi2_est_T3(i) = std( (x_obs-X*m_x) + sqrt(-1)*(y_obs-X*m_y) );
   final_variance_T3(i) = std( U1*m_x + sqrt(-1)*U1*m_y )/sqrt(2);
end

[val, index] = min(rms_error_T1);
fprintf('S=%d, T=1, stride=%d, min-u=%g, rms_error=%g, chi2-est=%g, u-variance=%g\n', S, stride, u_vals(index), val, (chi2_est_T1(index)/(position_error*sqrt(2)))^2, final_variance_T1(index) );
[val, index] = min(rms_error_T2);
fprintf('S=%d, T=2, stride=%d, min-a=%g, rms_error=%g, chi2-est=%g, a-variance=%g\n', S, stride, a_vals(index), val, (chi2_est_T2(index)/(position_error*sqrt(2)))^2, final_variance_T2(index) );
[val, index] = min(rms_error_T3);
fprintf('S=%d, T=3, stride=%d, min-j=%g, rms_error=%g, chi2-est=%g, j-variance=%g\n', S, stride, j_vals(index), val, (chi2_est_T3(index)/(position_error*sqrt(2)))^2, final_variance_T3(index) );

figure
scatter(u_vals, rms_error_T1);
xlog, ylog
vlines(u_true,'g--')

figure
scatter(a_vals, rms_error_T2);
xlog, ylog
vlines(a_true,'g--')

figure
scatter(j_vals, rms_error_T3);
xlog, ylog
vlines(j_true,'g--')


% Note that the rms errors reported below are for position, and not divided
% by the necessary sqrt(2).
%
% Measured std-u=0.195264, measured std-a=0.00210827, measured std-j=0.000421435
% Measured std-u=0.195264, measured std-a=2.1e-3, measured std-j=4.2e-4
% S=2, T=1, stride=1, length/10, min-u=0.19, rms-error: 2.946
% S=2, T=2, stride=1, length/10, min-a=0.003162, rms-error: 2.726

% Maximal tension parameters: u: 0.11 (rms: 4.12), a: 8e-4 (rms 3.48), j: 7.5e-6 (rms 3.47)
% S=3, T=1, stride=1, length/5, min-u=0.19, rms-error: 2.87
% S=3, T=2, stride=1, length/5, min-a=3.2e-3, rms-error: 2.56
% S=3, T=2, stride=1, length/5, min-j=3.7e-5, rms-error: 2.60

% S=2, T=1, stride=10, length, min-u=0.12, rms-error: 6.4, deduced-u=0.1956
% S=2, T=2, stride=10, length, min-a=7e-4, rms-error: 5.9, deduced-a=5e-4
% S=3, T=3, stride=10, length, min-j=5.2e-6, rms-error: 6.0, deduced-j=not possible

% S=2, T=1, stride=100, length, min-u=0.07, rms-error: 13.3, deduced-u=0.19
% S=2, T=2, stride=100, length, min-a=2e-4, rms-error: 12.6, deduced-a=2e-4
% S=3, T=3, stride=100, length, min-j=1.9e-6, rms-error: 12.8, deduced-a=6.3e-7

% Major conclusions:
% 1. The order of the spline doesn't seem to make a different (tried higher
% values than what is shown here.)
% 2. Tension definitely should be reduced as sampling becomes more sparse!
% 3. The order of the spline doesn't seem to matter much. Using
% accelerations was mildly better in this case, but not by much.
%
% I would argue that you use as high of order of tension as you can, before
% noise takes over. For well sampled data, this means 

% Looking at the ratio of 0.19^2/min-u^2 for the three subsampled cases:
% min-u decreases by 1.0, 2.5, 7.3 going from stride, 1, 10, 100
% min-a decreases by 1, 20, 250
% min-j decreases by 1, 50, 380

% No, I should do this from the measured std
% min-u decreases by 1.0, 2.5, 7.3 going from stride, 1, 10, 100
% min-a decreases by 0.4, 9, 110
% min-j decreases by 129, 6570, 50000

% Now using a better error metric that will account for wild variations
%
% (Note this has DF=10, so fewer knot points)
% Measured std-u=0.195264, measured std-a=0.00210827, measured std-j=0.000421435
% Gamma_u = 10
% S=3, T=1, stride=1, min-u=0.227585, rms_error=2.70874, chi2-est=0.930674, u-variance=0.187477
% S=3, T=2, stride=1, min-a=0.0026827, rms_error=2.47764, chi2-est=0.955499, a-variance=0.000472391
% S=3, T=3, stride=1, min-j=6.1054e-05, rms_error=2.50329, chi2-est=0.953498, j-variance=6.37223e-06
%
% Gamma_u = 1
% S=3, T=1, stride=10, min-u=0.108571, rms_error=6.45238, chi2-est=0.709578, u-variance=0.19165
% S=3, T=2, stride=10, min-a=0.001, rms_error=5.89656, chi2-est=0.809697, a-variance=0.000376931
% S=3, T=3, stride=10, min-j=8.48343e-06, rms_error=5.98009, chi2-est=0.850188, j-variance=2.17002e-06
% 
% Gamma_u = 0.1
% S=3, T=1, stride=100, min-u=0.0848343, rms_error=14.7327, chi2-est=0.121644, u-variance=0.187512
% S=3, T=2, stride=100, min-a=0.00026827, rms_error=14.2462, chi2-est=0.263045, a-variance=0.000238325
% S=3, T=3, stride=100, min-j=1.9307e-06, rms_error=14.3115, chi2-est=0.154203, j-variance=1.13312e-06
%
% min-u decreases by 1.0, 3, 5 going from stride, 1, 10, 100
% min-a decreases by 1.0, 7, 100
% min-j decreases by 1.0, 19, 380

% What are my big take homes for the day?
% 1. The deduced a and j variance from fits are much less than actual,
% presumably because we can't resolve the higher frequencies.
%
% 2. I cannot account for the required drop in tension for decreasing
% resolution in the signal. Actually, with u, maybe we can simply by
% looking at chi-squared.
%
% 3. The most consistent value appears to be the chi2 metric.
