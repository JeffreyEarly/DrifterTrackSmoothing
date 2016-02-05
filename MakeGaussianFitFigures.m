% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

drifters = load('sample_data/projected_ungridded_rho1_drifters.mat');

 distribution = 'student-t';
 distribution = 'gaussian';
 
 iDrifter = 6;

 x = drifters.x{iDrifter};
y = drifters.y{iDrifter};
t = drifters.t{iDrifter};

S = 3; % order of the spline
K = S+1;

maxlag = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Large error, large tension case
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma = 200; % error in meters
a = 2.5e-5;

sigma = 195;
a = 1.5e-5;

% optimal total error (and position error), given a
sigma = 195;
a = 1.5e-5;

% optimal acceleration error, given a
sigma = 65;
a = 1.5e-5;

position_pdf_big = @(z) exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi));
w = @(z)(sigma*sigma);

velocity_pdf_big = @(z) exp(-(z.*z)/(2*a*a))/(a*sqrt(2*pi));
gamma = [0; 1/a^2; 0];
    
dx = ones(size(x))*sigma;
dy = ones(size(y))*sigma;

[m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = bspline_bivariate_fit_with_tension(t,x,y,dx,dy,S,gamma, w);

Xq = squeeze(Bq(:,:,1));
x_fit_big = Xq*m_x;
y_fit_big = Xq*m_y;

a_big = [squeeze(Bq(:,:,3))*m_x; squeeze(Bq(:,:,3))*m_y];

X = squeeze(B(:,:,1));
error_x_big = X*m_x - x;
error_y_big = X*m_y - y;
error_big = [error_x_big;error_y_big];
ACx_big = Autocorrelation(error_x_big,maxlag);
ACy_big = Autocorrelation(error_y_big,maxlag);
AC_big = (ACx_big + ACy_big)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Small error, small tension case
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma = 20; % error in meters
a = 9e-5;

% optimal total error, both parameters free
sigma = 80;
a = 1.1826e-04;

% optimal total error, fixed sigma
sigma = 200;
a = 1.8125e-6;

% optimal acceleration error, both parameters free
sigma = 131.8940;
a = 10^(-5.6154);

position_pdf_small = @(z) exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi));
w = @(z)(sigma*sigma);
velocity_pdf_small = @(z) exp(-(z.*z)/(2*a*a))/(a*sqrt(2*pi));
gamma = [0; 1/a^2; 0];
    
dx = ones(size(x))*sigma;
dy = ones(size(y))*sigma;

[m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = bspline_bivariate_fit_with_tension(t,x,y,dx,dy,S,gamma, w);

Xq = squeeze(Bq(:,:,1));
x_fit_small = Xq*m_x;
y_fit_small = Xq*m_y;

a_small = [squeeze(Bq(:,:,3))*m_x; squeeze(Bq(:,:,3))*m_y];

std(a_small)

X = squeeze(B(:,:,1));
error_x_small = X*m_x - x;
error_y_small = X*m_y - y;
error_small = [error_x_small;error_y_small];
ACx_small = Autocorrelation(error_x_small,maxlag);
ACy_small = Autocorrelation(error_y_small,maxlag);
AC_small = (ACx_small + ACy_small)/2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Position fit figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(2,2,[1 3])
s = 1/1000;
plot(s*x_fit_small,s*y_fit_small,'b'), hold on
plot(s*x_fit_big,s*y_fit_big,'k')
scatter(s*x,s*y,5)
xlabel('x (km)')
ylabel('y (km)')

subplot(2,2,2)
plot(tq/3600,s*x_fit_small,'b'), hold on
plot(tq/3600,s*x_fit_big,'k')
scatter(drifters.t{iDrifter}/3600,s*x,5)
xlabel('t (hours)')
ylabel('x (km)')

subplot(2,2,4)
plot(tq/3600,s*y_fit_small,'b'), hold on
plot(tq/3600,s*y_fit_big,'k')
scatter(drifters.t{iDrifter}/3600,s*y,5)
xlabel('t (hours)')
ylabel('y (km)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Position error histogram
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

subplot(2,2,1)
plot_hist_with_pdf( error_big, position_pdf_big, 1000, 50 )

subplot(2,2,2)
plot_hist_with_pdf( error_small, position_pdf_small, 100, 50 )

subplot(2,2,3)
plot_hist_with_pdf( a_big, velocity_pdf_big, 1e-4, 50 )

subplot(2,2,4)
plot_hist_with_pdf( a_small, velocity_pdf_small, 5e-4, 50 )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Autocorrelation sequence
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

plot(AC_small), hold on
plot(AC_big)



