% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

drifters = load('sample_data/projected_ungridded_rho1_drifters.mat');

% Drifter to highlight in the final plots
choiceDrifter = 6;
Ndrifters = length(drifters.x);

S = 3; % order of the spline
K = S+1;

maxlag = 30;

% How many data points do we have total
Nall = 0;
for iDrifter = 1:9
    Nall = Nall + length(drifters.x{iDrifter});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Large error, large tension case
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nu = 2.00; sigma = 8;
a = 0.8e-5;

position_pdf_big = @(z) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2));
w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));

velocity_pdf_big = @(z) exp(-(z.*z)/(2*a*a))/(a*sqrt(2*pi));
velocity_cdf_big = @(z) 0.5*(1 + erf(z/(a*sqrt(2))));    

ACx_big = zeros(maxlag+1,1);
ACy_big = zeros(maxlag+1,1);
a_big = [];
error_big = zeros(2*Nall,1);
iEpsilon = 1;
for iDrifter = 1:Ndrifters
    x = drifters.x{iDrifter};
    y = drifters.y{iDrifter};
    t = drifters.t{iDrifter};
    N = length(t);
    
    dx = ones(size(x))*sigma;
    dy = ones(size(y))*sigma;
    [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = bspline_bivariate_fit_with_tension(t,x,y,dx,dy,S,[0; 1/a^2; 0], w);
    
    Xq = squeeze(Bq(:,:,1));
    if (iDrifter == choiceDrifter)
        x_fit_big = Xq*m_x;
        y_fit_big = Xq*m_y;
    end
    
    a_big = [a_big; squeeze(Bq(:,:,3))*m_x; squeeze(Bq(:,:,3))*m_y];
    
    X = squeeze(B(:,:,1));
    error_x_big = X*m_x - x;
    error_y_big = X*m_y - y;
    error_big(iEpsilon:(iEpsilon+2*N-1)) = [error_x_big;error_y_big];
    iEpsilon = iEpsilon + 2*N + 1;
    
    ACx_big = ACx_big + Autocorrelation(error_x_big,maxlag);
    ACy_big = ACy_big + Autocorrelation(error_y_big,maxlag);
end

ACx_big = ACx_big/Ndrifters;
ACy_big = ACy_big/Ndrifters;
AC_big = (ACx_big + ACy_big)/2;

log10(std(a_big)/a)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Small error, small tension case
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Very good set of parameters
nu = 2.0; sigma = 4;
a = 1.6e-5;

% Very good set of parameters
nu = 2.0; sigma = 20;
a = 0.52e-5;


position_pdf_small = @(z) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2));
w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));
velocity_pdf_small = @(z) exp(-(z.*z)/(2*a*a))/(a*sqrt(2*pi));
velocity_cdf_small = @(z) 0.5*(1 + erf(z/(a*sqrt(2))));    

ACx_small = zeros(maxlag+1,1);
ACy_small = zeros(maxlag+1,1);
a_small = [];
error_small = zeros(2*Nall,1);
iEpsilon = 1;
for iDrifter = 1:Ndrifters
    x = drifters.x{iDrifter};
    y = drifters.y{iDrifter};
    t = drifters.t{iDrifter};
    N = length(t);
    
    dx = ones(size(x))*sigma;
    dy = ones(size(y))*sigma;
    [m_x,m_y,Cm_x,Cm_y,B,Bq,tq2] = bspline_bivariate_fit_with_tension(t,x,y,dx,dy,S,[0; 1/a^2; 0], w);
    
    Xq = squeeze(Bq(:,:,1));
    if (iDrifter == choiceDrifter)
        tq = tq2;
        x_fit_small = Xq*m_x;
        y_fit_small = Xq*m_y;
    end
    
    a_small = [a_small; squeeze(Bq(:,:,3))*m_x; squeeze(Bq(:,:,3))*m_y];
    
    X = squeeze(B(:,:,1));
    error_x_small = X*m_x - x;
    error_y_small = X*m_y - y;
    error_small(iEpsilon:(iEpsilon+2*N-1)) = [error_x_small;error_y_small];
    iEpsilon = iEpsilon + 2*N + 1;
    
    ACx_small = ACx_small + Autocorrelation(error_x_small,maxlag);
    ACy_small = ACy_small + Autocorrelation(error_y_small,maxlag);
end

ACx_small = ACx_small/Ndrifters;
ACy_small = ACy_small/Ndrifters;
AC_small = (ACx_small + ACy_small)/2;

log10(std(a_small)/a)

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
scatter(s*drifters.x{choiceDrifter},s*drifters.y{choiceDrifter},5)
xlabel('x (km)')
ylabel('y (km)')

subplot(2,2,2)
plot(tq/3600,s*x_fit_small,'b'), hold on
plot(tq/3600,s*x_fit_big,'k')
scatter(drifters.t{choiceDrifter}/3600,s*drifters.x{choiceDrifter},5)
xlabel('t (hours)')
ylabel('x (km)')

subplot(2,2,4)
plot(tq/3600,s*y_fit_small,'b'), hold on
plot(tq/3600,s*y_fit_big,'k')
scatter(drifters.t{choiceDrifter}/3600,s*drifters.y{choiceDrifter},5)
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
plot_hist_with_pdf( error_small, position_pdf_small, 100, 50 )

subplot(2,2,3)
plot_hist_with_pdf( a_big, velocity_pdf_big, 5e-5, 50 )

subplot(2,2,4)
plot_hist_with_pdf( a_small, velocity_pdf_small, 5e-5, 50 )

figure
subplot(1,2,1)
plot_hist_with_cdf( a_big, velocity_cdf_big, 5e-5, 50 )

subplot(1,2,2)
plot_hist_with_cdf( a_small, velocity_cdf_small, 10e-5, 50 )

x = sort(a_small-mean(a_small));
n = length(x);
y_data = (1:n)'/n;
y = velocity_cdf_small(x);
figure
plot(x,y_data), hold on, plot(x,y)
D = max(abs(y-y_data))

lambda = (sqrt(n) + 0.12 + 0.11/sqrt(n))*D;
j=1:25;
p = sum(2*((-1).^(j-1)).*exp(-2*j.*j*lambda*lambda))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Autocorrelation sequence
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

plot(AC_small), hold on
plot(AC_big)



