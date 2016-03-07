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

stride = 1;
indices = 1:stride:floor(length(t)/5);
x_obs = x(indices) + epsilon_x(indices);
y_obs = y(indices) + epsilon_y(indices);
t_obs = t(indices);
sigma = position_error;
S = 3;
T = 1;
K = S+1;

% t_knot = NaturalKnotsForSpline( t, K, 1 );
% B_all = bspline(t,t_knot,K);

u_vals = 10.^(linspace(-2,1,15))';
rms_error_T1 = zeros(size(u_vals));
for i = 1:length(u_vals)
   [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t_obs,x_obs,y_obs,sigma,S,T,u_vals(i));
   X = squeeze(B(:,:,1));
   rms_error_T1(i) = std( (x(indices)-X*m_x) + sqrt(-1)*(y(indices)-X*m_y) );   
end

T = 2;
a_vals = 10.^(linspace(-4,-1,15))';
rms_error_T2 = zeros(size(a_vals));
for i = 1:length(a_vals)
   [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t_obs,x_obs,y_obs,sigma,S,T,a_vals(i));
   X = squeeze(B(:,:,1));
   rms_error_T2(i) = std( (x(indices)-X*m_x) + sqrt(-1)*(y(indices)-X*m_y) );   
end

T = 3;
j_vals = 10.^(linspace(-7,-3,15))';
rms_error_T3 = zeros(size(j_vals));
for i = 1:length(j_vals)
   [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t_obs,x_obs,y_obs,sigma,S,T,j_vals(i));
   X = squeeze(B(:,:,1));
   rms_error_T3(i) = std( (x(indices)-X*m_x) + sqrt(-1)*(y(indices)-X*m_y) );   
end

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
