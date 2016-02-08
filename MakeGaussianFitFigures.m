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
sigma = 200; a = 2.5e-5;

sigma = 195; a = 1.5e-5;

% optimal total error (and position error), given a
sigma = 195; a = 1.5e-5;

% optimal acceleration error, given a
sigma = 65; a = 1.5e-5;
sigma = 65; a = 0.5e-5;

position_pdf_big = @(z) exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi));
w = @(z)(sigma*sigma);

velocity_pdf_big = @(z) exp(-(z.*z)/(2*a*a))/(a*sqrt(2*pi));
velocity_cdf_big = @(z) 0.5*(1 + erf(z/(a*sqrt(2))));    

ACx_big = zeros(maxlag+1,1);
ACy_big = zeros(maxlag+1,1);
a_big = [];
ax_big = [];
ay_big = [];
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
    ax_big = [ax_big; squeeze(Bq(:,:,3))*m_x];
    ay_big = [ay_big; squeeze(Bq(:,:,3))*m_y];
    
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Small error, small tension case
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma = 20; % error in meters
a = 9e-5;

% optimal total error, both parameters free
% sigma = 80;
% a = 1.1826e-04;
% 
% % optimal total error, fixed sigma
% sigma = 200;
% a = 1.8125e-6;
% 
% % optimal acceleration error, both parameters free
% sigma = 131.8940;
% a = 10^(-5.6154);
% a = 1e-3;

position_pdf_small = @(z) exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi));
w = @(z)(sigma*sigma);
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Position fit figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

subplot(2,2,[1 3])
s = 1/1000;
plot(s*x_fit_small,s*y_fit_small, 'LineWidth', 0.5*scaleFactor, 'Color', 'b'), hold on
plot(s*x_fit_big,s*y_fit_big, 'LineWidth', 1.0*scaleFactor, 'Color','k')
scatter(s*drifters.x{choiceDrifter},s*drifters.y{choiceDrifter},(5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('x (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('y (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
% xlim([-10 -7.4])
% ylim([13.5 19])

sp1 = subplot(2,2,2);
plot(tq/3600,s*x_fit_small, 'LineWidth', 0.5*scaleFactor, 'Color','b'), hold on
plot(tq/3600,s*x_fit_big, 'LineWidth', 1.0*scaleFactor, 'Color','k')
scatter(drifters.t{choiceDrifter}/3600,s*drifters.x{choiceDrifter},(5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
ylabel('x (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
set(sp1,'yaxislocation','right')
% xlim([min(t) 40])
% ylim([-10 -7.4])

sp1 = subplot(2,2,4)
plot(tq/3600,s*y_fit_small, 'LineWidth', 0.5*scaleFactor, 'Color','b'), hold on
plot(tq/3600,s*y_fit_big, 'LineWidth', 1.0*scaleFactor, 'Color','k')
scatter(drifters.t{choiceDrifter}/3600,s*drifters.y{choiceDrifter},(5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('t (hours)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('y (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
set(sp1,'yaxislocation','right')
% xlim([min(t) 40])
% ylim([13.5 19])


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
plot(tq/3600,s*x_fit_big, 'LineWidth', 0.5*scaleFactor, 'Color',0.0*[1.0 1.0 1.0])
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

print('-depsc2', 'figures/gaussianfit.eps')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Autocorrelation sequence
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FigureSize = [50 50 figure_width_1col+6 135*scaleFactor];

fig1 = figure('Units', 'points', 'Position', FigureSize)
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

plot(AC_small, 'LineWidth', 0.5*scaleFactor, 'Color',0.4*[1.0 1.0 1.0]), hold on
plot(AC_big, 'LineWidth', 0.5*scaleFactor, 'Color',0.0*[1.0 1.0 1.0])
xlabel('lag', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('autocorrelation', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylim([-1 1])
xlim([0 25])

fig1 = tightfig;
fig1.Position = FigureSize;
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];
fig1.PaperPositionMode = 'auto';

print('-depsc2', 'figures/gaussianacf.eps')

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
subplot(2,1,1)
plot(AC_small), hold on
plot(AC_big)
subplot(2,1,2)
[p1, Q1] = LjungBoxTest(AC_small, Nall);
[p2, Q2] = LjungBoxTest(AC_big, Nall);
plot(p1), hold on
plot(p2)




