% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

drifters = load('sample_data/rho1_drifters_projected_ungridded.mat');

% Drifter to highlight in the final plots
choiceDrifter = 6;

S = 3; % order of the spline
K = S+1;
T = 2; % order of the tension

% Pull out the data of interest
x = drifters.x{choiceDrifter};
y = drifters.y{choiceDrifter};
t = drifters.t{choiceDrifter};
N = length(t);

% Standard deviation of the noise
sigma = 10;

% Estimate the velocities...
u_estimate_spectral = EstimateRMSVelocityFromSpectrum(t,x,sigma, 1);
v_estimate_spectral = EstimateRMSVelocityFromSpectrum(t,y,sigma, 1);

% ...and accelerations
ax_estimate_spectral = EstimateRMSAccelerationFromSpectrum(t,x,sigma);
ay_estimate_spectral = EstimateRMSAccelerationFromSpectrum(t,y,sigma);

% Now estimate the optimal tension parameter
dt = t(2)-t(1);
expectedDOFx = 1 + 3*sigma/(u_estimate_spectral*dt);
lambda_x = (expectedDOFx-1)/(expectedDOFx*ax_estimate_spectral^2);

expectedDOFy = 1 + 3*sigma/(u_estimate_spectral*dt);
lambda_y = (expectedDOFy-1)/(expectedDOFy*ay_estimate_spectral^2);

DF = 1;
[m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t,x,y,sigma,S,T,lambda_x,lambda_y,DF);
% X = squeeze(B(:,:,1));

Xq = squeeze(Bq(:,:,1));
x_fit_big = Xq*m_x;
y_fit_big = Xq*m_y;

[m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t,x,y,sigma,S,T,1e7*lambda_x,1e7*lambda_y,DF);

Xq = squeeze(Bq(:,:,1));
x_fit_small = Xq*m_x;
y_fit_small = Xq*m_y;

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

sp1 = subplot(2,2,4);
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

fig1 = figure('Units', 'points', 'Position', FigureSize);
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
ylim([5.58 9.18])

% packfig(2,2)
fig1 = tightfig;
fig1.Position = FigureSize;
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];
fig1.PaperPositionMode = 'auto';

return

% print('-depsc2', 'figures/gaussianfit.eps')

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

plot( (1:length(AC_small))'-1, AC_small, 'LineWidth', 0.5*scaleFactor, 'Color',0.4*[1.0 1.0 1.0]), hold on
plot( (1:length(AC_big))'-1, AC_big, 'LineWidth', 0.5*scaleFactor, 'Color',0.0*[1.0 1.0 1.0])
xlabel('lag', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('autocorrelation', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylim([-1 1])
xlim([0 15])

fig1 = tightfig;
fig1.Position = FigureSize;
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];
fig1.PaperPositionMode = 'auto';

% print('-depsc2', 'figures/gaussianacf.eps')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Position error histogram
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

subplot(2,2,1)
plot_hist_with_pdf( error_big, position_pdf_big, 60, 100 )

subplot(2,2,2)
plot_hist_with_pdf( error_small, position_pdf_small, 60, 100 )

subplot(2,2,3)
plot_hist_with_pdf( a_big, velocity_pdf_big, 1e-4, 50 )

subplot(2,2,4)
plot_hist_with_pdf( a_small, velocity_pdf_small, 5e-4, 50 )

figure
plot_hist_with_pdf( error_big, position_pdf_big, 200, 50, 0.3*[1.0 1.0 1.0] )
plot_hist_with_pdf( error_small, position_pdf_small, 200, 50, 0.3*[1.0 1.0 1.0] )

% figure
% plot_hist_with_pdf( a_big, velocity_pdf_big, 1e-4, 50, 0.3*[1.0 1.0 1.0] )
% plot_hist_with_pdf( a_small, velocity_pdf_small, 1e-4, 50, 0.3*[1.0 1.0 1.0] )

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




