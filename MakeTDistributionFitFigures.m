% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

drifters = load('sample_data/projected_ungridded_rho1_drifters.mat');

% Drifter to highlight in the final plots
choiceDrifter = 6;
Ndrifters = length(drifters.x);

S = 3; % order of the spline
K = S+1;
sigma_gps = 1.75;

shouldDiscardOutliersInACF = 0;
shouldUseSignedACF = 0;

maxlag = 30;

% How many data points do we have total
Nall = 0;
for iDrifter = 1:Ndrifters
    Nall = Nall + length(drifters.x{iDrifter});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Large error, large tension case
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = 3; % order of the spline
K = S+1;
T = 2; % order at which tension is applied

nu = 2.00; sigma = 8;
a = 0.8e-5; % input acceleration variances matches output
a= 1.1994e-05; % Ljung-Box minimum

nu = 2.00; sigma = 10.0; a = 1.2317e-05;
nu = 2.00; sigma = 5.0; a = 1.1272e-05;
nu = 2.00; sigma = 1.25; a = 7.2975e-06;
nu = 2.00; sigma = 1.25; a = 6.9e-06;

nu = 2.00; sigma = 4; a = 9.2956e-06; % sigma_gps = 1.75 min in the KS test
% nu = 2.00; sigma = 4; a = 5.1614e-06; % sigma_gps = 4 min in the KS test
% nu = 2.00; sigma = 3*sigma_gps; a =1.29e-05; % sigma_gps = 1.75 min in the KS test

% S = 4; K = S+1;
% nu = 2.0; sigma = 4.0; T=2; a = 9e-06;
% 
% S = 5; K = S+1;
% nu = 2.0; sigma = 4.0; T=3; a = 3*1.3e-09;

% S = 5; K = S+1;
% nu = 2.0; sigma = 4.0; T=4; a = 2*1.3e-12;

% S = 6; K = S+1;
% nu = 2.0; sigma = 4.0; T=5; a = 1.4e-15;

S = 4; K = S+1;
nu = 10.0; sigma = sigma_gps; T=3; a = 5.9113e-09;

S = 10; K = S+1;
nu = 10.0; sigma = sigma_gps; T=2; a = 5.16e-06;

S = 3; K = S+1;
nu = 5.0; sigma = sigma_gps; T=2; a = 5e-06;

S = 3; K = S+1; sigma_gps = 6.0;
nu = 6.0; sigma = sigma_gps; T=2; a = 5.3462e-06;

gaussian_pdf_big = @(z) exp(-(z.*z)/(2*sigma_gps*sigma_gps))/(sigma_gps*sqrt(2*pi));

position_pdf_big = @(z) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2));
w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));

velocity_pdf_big = @(z) exp(-(z.*z)/(2*a*a))/(a*sqrt(2*pi));
velocity_cdf_big = @(z) 0.5*(1 + erf(z/(a*sqrt(2))));    

ACx_big = zeros(maxlag+1,1);
ACy_big = zeros(maxlag+1,1);
a_big = [];
ax_big = [];
ay_big = [];
n_acf_big = 0;
% error_big = zeros(2*Nall,1);
% dist_error_big  = zeros(Nall,1);
error_big = [];
dist_error_big = [];
iEpsilon = 1;
iDistError = 1;
for iDrifter = 1:Ndrifters
    x = drifters.x{iDrifter};
    y = drifters.y{iDrifter};
    t = drifters.t{iDrifter};
    N = length(t);
    
    dx = ones(size(x))*sigma;
    dy = ones(size(y))*sigma;
    tension = zeros(S,1);
    tension(T) = 1/a^2;
    [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = bspline_bivariate_fit_with_tension(t,x,y,dx,dy,S,tension, w);
    
    X = squeeze(B(:,:,1));
    error_x_big = X*m_x - x;
    error_y_big = X*m_y - y;
    dist_error = sqrt( error_x_big.*error_x_big + error_y_big.*error_y_big );
    badNuts = dist_error > 50;
    
    t(badNuts) = [];
    x(badNuts) = [];
    y(badNuts) = [];
    
    sprintf('Decreased points by %d\n',sum(badNuts))
    [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = bspline_bivariate_fit_with_tension(t,x,y,ones(size(x))*sigma,ones(size(x))*sigma,S, 2*tension, w);
    
    Xq = squeeze(Bq(:,:,1));
    if (iDrifter == choiceDrifter)
        x_fit_big = Xq*m_x;
        y_fit_big = Xq*m_y;
        tq_big = tq;
    end
    
    a_big = [a_big; squeeze(Bq(:,:,T+1))*m_x; squeeze(Bq(:,:,T+1))*m_y];
    ax_big = [ax_big; squeeze(Bq(:,:,T+1))*m_x];
    ay_big = [ay_big; squeeze(Bq(:,:,T+1))*m_y];
    
    X = squeeze(B(:,:,1));
    error_x_big = X*m_x - x;
    error_y_big = X*m_y - y;
%     error_big(iEpsilon:(iEpsilon+2*N-1)) = [error_x_big;error_y_big];
%     dist_error_big(iDistError:(iDistError+N-1)) = sqrt( error_x_big.*error_x_big + error_y_big.*error_y_big );
%     iEpsilon = iEpsilon + 2*N;
%     iDistError = iDistError + N;
    
    error_big = [error_big; error_x_big; error_y_big];
    dist_error_big = [dist_error_big; sqrt( error_x_big.*error_x_big + error_y_big.*error_y_big )];
    
    if shouldUseSignedACF == 1
        error_x_big = sign(error_x_big);
        error_y_big = sign(error_y_big);
    end
    
    if shouldDiscardOutliersInACF == 1
        error_x_big(abs(error_x_big) > 6*sigma) = [];
        error_y_big(abs(error_y_big) > 6*sigma) = [];
    end
    
    ACx_big = ACx_big + Autocorrelation(error_x_big,maxlag);
    ACy_big = ACy_big + Autocorrelation(error_y_big,maxlag);
    n_acf_big = n_acf_big + length(error_x_big) + length(error_y_big);
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

S = 3; % order of the spline
K = S+1;


% Very good set of parameters
nu = 2.0; sigma = 4;
a = 1.6e-5;
a = 1.0799e-05; % Ljung-Box test

% Very good set of parameters
% nu = 2.0; sigma = 20;
% a = 0.52e-5;

% Optimal for the Ljung-Box test
nu = 2.0; sigma = 1.7; a = 9.3254e-06;

% Nudge down the acceleration to reject outliers.
nu = 2.0; sigma = 1.7; a = 8e-06;

% Optimal for the Ljung-Box test
nu = 2.0; sigma = 2.5; a = 9e-06;

% Nudge down the acceleration to reject outliers. %% Current favorite!!!
S = 3; K = S+1;
nu = 2.0; sigma = 4.0; a = 7e-06;

S = 3; K = S+1;
nu = 2.00; sigma = 4; a = 9.2956e-06; % sigma_gps = 1.75 min in the KS test
nu = 2.00; sigma = 4; a = 5.1614e-06; % sigma_gps = 4 min in the KS test
nu = 2.00; sigma = 3*sigma_gps; a =1.29e-05; % sigma_gps = 1.75 min in the KS test
nu = 5; sigma=sigma_gps; a = 4.1839e-06;

sigma_gps = 1.0;
nu = 6; sigma=sigma_gps; a = 2.2412e-06;

S = 3; K = S+1; sigma_gps = 6.0;
nu = 6.0; sigma = sigma_gps; T=2; a = 5.3462e-06;

% nu = 10; sigma=sigma_gps; a = 4.7653e-06;


% Super high order spline, lets us lower the tolerance on a
% S = 20; K = S+1;
% nu = 2.0; sigma = 4.0; a = 9e-06;

% ACF experiment
% nu = 2.0; sigma = 4.0; a = 3e-06;



gaussian_pdf_small = @(z) exp(-(z.*z)/(2*sigma_gps*sigma_gps))/(sigma_gps*sqrt(2*pi));

position_pdf_small = @(z) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2));
w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));

velocity_pdf_small = @(z) exp(-(z.*z)/(2*a*a))/(a*sqrt(2*pi));
exponential_pdf_small = @(z) exp(-(abs(z))/(a))/(2*a);
velocity_cdf_small = @(z) 0.5*(1 + erf(z/(a*sqrt(2))));    

ACx_small = zeros(maxlag+1,1);
ACy_small = zeros(maxlag+1,1);
a_small = [];
ax_small = [];
ay_small = [];
n_acf_small = 0;
error_small = zeros(2*Nall,1);
iEpsilon = 1;
for iDrifter = 1:Ndrifters
    x = drifters.x{iDrifter};
    y = drifters.y{iDrifter};
    t = drifters.t{iDrifter};
    N = length(t);
    
    dx = ones(size(x))*sigma;
    dy = ones(size(y))*sigma;
    tension = zeros(S,1);
    tension(2) = 1/a^2;
    [m_x,m_y,Cm_x,Cm_y,B,Bq,tq2] = bspline_bivariate_fit_with_tension(t,x,y,dx,dy,S,tension, w);
    
    Xq = squeeze(Bq(:,:,1));
    if (iDrifter == choiceDrifter)
        tq = tq2;
        x_fit_small = Xq*m_x;
        y_fit_small = Xq*m_y;
    end
    
    a_small = [a_small; squeeze(Bq(:,:,3))*m_x; squeeze(Bq(:,:,3))*m_y];
    ax_small = [ax_small; squeeze(Bq(:,:,3))*m_x];
    ay_small = [ay_small; squeeze(Bq(:,:,3))*m_y];
    
    X = squeeze(B(:,:,1));
    error_x_small = X*m_x - x;
    error_y_small = X*m_y - y;
    error_small(iEpsilon:(iEpsilon+2*N-1)) = [error_x_small;error_y_small];
    iEpsilon = iEpsilon + 2*N;
    
    if shouldUseSignedACF == 1
        error_x_small = sign(error_x_small);
        error_y_small = sign(error_y_small);
    end
        
    if shouldDiscardOutliersInACF == 1
        error_x_small(abs(error_x_small) > 6*sigma) = [];
        error_y_small(abs(error_y_small) > 6*sigma) = [];
    end
    
    n_acf_small = n_acf_small + length(error_x_small) + length(error_y_small);
    
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
plot(tq_big/3600,s*x_fit_big,'k')
scatter(drifters.t{choiceDrifter}/3600,s*drifters.x{choiceDrifter},5)
xlabel('t (hours)')
ylabel('x (km)')

subplot(2,2,4)
plot(tq/3600,s*y_fit_small,'b'), hold on
plot(tq_big/3600,s*y_fit_big,'k')
scatter(drifters.t{choiceDrifter}/3600,s*drifters.y{choiceDrifter},5)
xlabel('t (hours)')
ylabel('y (km)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Position fit figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FigureSize = [50 50 figure_width_2col+8 150*scaleFactor];

fig1 = figure('Units', 'points', 'Position', FigureSize)
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

plot(tq/3600,s*x_fit_small, 'LineWidth', 0.5*scaleFactor, 'Color',0.4*[1.0 1.0 1.0]), hold on
plot(tq_big/3600,s*x_fit_big, 'LineWidth', 0.5*scaleFactor, 'Color',0.0*[1.0 1.0 1.0])
scatter(drifters.t{choiceDrifter}/3600,s*drifters.x{choiceDrifter},(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('t (hours)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('x (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)

xlim([124 149])
ylim([4.4 8.0])

% packfig(2,2)
fig1 = tightfig;
fig1.Position = FigureSize;
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];
fig1.PaperPositionMode = 'auto';

print('-depsc2', 'figures/tdistributionfit.eps')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Position error histogram
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

subplot(2,2,1)
plot_hist_with_pdf( error_big, position_pdf_big, 60, 100 )
% plot_hist_with_pdf( error_big, gaussian_pdf_big, 60, 100 )

subplot(2,2,2)
plot_hist_with_pdf( error_small, position_pdf_small, 60, 100 )
% plot_hist_with_pdf( error_small, gaussian_pdf_small, 20, 100 )

subplot(2,2,3)
plot_hist_with_pdf( a_big, velocity_pdf_big, 10e-5, 50 )

subplot(2,2,4)
plot_hist_with_pdf( a_small, velocity_pdf_small, 10e-5, 50 )
% plot_hist_with_pdf( a_small, exponential_pdf_small, 10e-5, 50 )


figure
subplot(1,2,1)
plot_hist_with_cdf( a_big, velocity_cdf_big, 10e-5, 50 )

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
subplot(2,1,1)
plot(AC_small), hold on
plot(AC_big)
subplot(2,1,2)
[p1, Q1] = LjungBoxTest(AC_small, n_acf_small);
[p2, Q2] = LjungBoxTest(AC_big, n_acf_big);
plot(p1), hold on
plot(p2)


