% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

drifters = load('sample_data/projected_ungridded_rho1_drifters.mat');

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
nu = 2.03; sigma = 20;
a = 2.5e-5;


nu = 2.0; sigma = 3;
a = 1e-5;

nu = 2.005; sigma = 8;
a = 0.5e-5;

position_pdf_big = @(z) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2));
w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));
velocity_pdf_big = @(z) exp(-(z.*z)/(2*a*a))/(a*sqrt(2*pi));
    
dx = ones(size(x))*sigma;
dy = ones(size(y))*sigma;

[m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = bspline_bivariate_fit_with_tension(t,x,y,dx,dy,S,[0; 1/a^2; 0], w);

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
nu = 2.005; sigma = 8;
position_pdf_small = @(z) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2));
w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));
a = 1.5e-5;
velocity_pdf_small = @(z) exp(-(z.*z)/(2*a*a))/(a*sqrt(2*pi));
    
dx = ones(size(x))*sigma;
dy = ones(size(y))*sigma;

[m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = bspline_bivariate_fit_with_tension(t,x,y,dx,dy,S,[0; 1/a^2; 0], w);

Xq = squeeze(Bq(:,:,1));
x_fit_small = Xq*m_x;
y_fit_small = Xq*m_y;

a_small = [squeeze(Bq(:,:,3))*m_x; squeeze(Bq(:,:,3))*m_y];

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
plot_hist_with_pdf( error_big, position_pdf_big, 50, 50 )

subplot(2,2,2)
plot_hist_with_pdf( error_small, position_pdf_small, 50, 50 )

subplot(2,2,3)
plot_hist_with_pdf( a_big, velocity_pdf_big, 1e-4, 50 )

subplot(2,2,4)
plot_hist_with_pdf( a_small, velocity_pdf_small, 1e-4, 50 )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Autocorrelation sequence
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

plot(AC_small), hold on
plot(AC_big)



