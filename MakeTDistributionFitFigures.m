% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults
addpath('support')

% Drifter to highlight in the final plots
choiceDrifter = 6;

shouldSaveFigures = 1;

% How many data points do we have total
drifters = load('smoothed_interpolated_rho1_drifters_T2');
load('smoothed_interpolated_rho1_drifters_T2');
sigma_t = sigma;
Ndrifters = length(drifters.x);

t_pdf = @(z) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2));

velocity_pdf = @(z) exp(-(z.*z)/(2*a*a))/(a*sqrt(2*pi));
velocity_cdf = @(z) 0.5*(1 + erf(z/(a*sqrt(2))));    

a_big = [];
ax_big = [];
ay_big = [];
error = [];
dist_error = [];
error_despiked = [];
dist_error_despiked = [];
for iDrifter = 1:Ndrifters
    x = drifters.x{iDrifter};
    y = drifters.y{iDrifter};
    t = drifters.t{iDrifter};
    N = length(t);
  
    a_big = [a_big; drifters.ax{iDrifter};  drifters.ay{iDrifter}];
    ax_big = [ax_big; drifters.ax{iDrifter}];
    ay_big = [ay_big; drifters.ay{iDrifter}];
    
    x_error = drifters.x_error{iDrifter};
    y_error = drifters.y_error{iDrifter};
    
    error = [error; x_error; y_error];
    dist = sqrt( x_error.*x_error + y_error.*y_error );
    dist_error = [dist_error; dist];
    
    if iDrifter == choiceDrifter
        rejectedPointIndices = find(dist > outlierCut);
    end
    
    x_error = drifters.x_error_despiked{iDrifter};
    y_error = drifters.y_error_despiked{iDrifter};
    
    error_despiked = [error_despiked; x_error; y_error];
    dist_error_despiked = [dist_error_despiked; sqrt( x_error.*x_error + y_error.*y_error )];
end



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

s = 1/1000; % scale
plot(drifters.t{choiceDrifter}/3600,s*drifters.x{choiceDrifter}, 'LineWidth', 0.5*scaleFactor, 'Color',0.4*[1.0 1.0 1.0]), hold on
scatter(drifters.t_raw{choiceDrifter}(rejectedPointIndices)/3600,s*drifters.x_raw{choiceDrifter}(rejectedPointIndices),(6.5*scaleFactor)^2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
scatter(drifters.t_raw{choiceDrifter}/3600,s*drifters.x_raw{choiceDrifter},(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
xlabel('t (hours)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('x (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)

xlim([124 149])
ylim([5.58 9.18])

fig1 = tightfig;
fig1.Position = FigureSize;
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];
fig1.PaperPositionMode = 'auto';

if shouldSaveFigures == 1
print('-depsc2', 'figures/tdistributionfit.eps')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Error distribution figures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FigureSize = [50 50 figure_width_1col+4 150*scaleFactor];

fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error PDF plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = error_despiked;
histwidth = 80;
nbins = 100;

% Create the bins for the data
binwidth = histwidth/nbins;
edges = [-histwidth*100;((-histwidth/2+binwidth):binwidth:(histwidth/2-binwidth))';histwidth*100];
binleft = linspace((-histwidth/2),(histwidth/2-binwidth),nbins)';

% plot the data

% This is the data that includes the outliers
if 0
    countRaw = histcounts(error,edges)';
    % countRaw(2:end-1) = 0;
    g = bar(binleft, countRaw/(length(error)*binwidth), 'histc'); hold on;
    g.FaceColor = 0.2*[1.0 1.0 1.0];
end

% this is the data that doesn't
count = histcounts(data,edges)';
g = bar(binleft, count/(length(data)*binwidth), 'histc'); hold on;
g.FaceColor = 0.8*[1.0 1.0 1.0];

% create bins for the analytical pdf
xi_left = linspace(-histwidth/2,-histwidth/2+binwidth,10)';
xi_mid = linspace(-histwidth/2+binwidth,histwidth/2-binwidth,100)';
xi_right = linspace(histwidth/2-binwidth,histwidth/2,10)';
xi = [xi_left;xi_mid;xi_right];

% plot the analytical pdf
pdf = t_pdf;
edgedensity = integral(pdf,(histwidth/2-binwidth),2*histwidth)/binwidth;
denfunc = edgedensity*ones(size(xi));
denfunc(11:110) = pdf(xi_mid);
plot(xi,denfunc,'LineWidth',2*scaleFactor,'Color',0.0*[1.0 1.0 1.0])

xlabel('meters', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
set( gca, 'FontSize', figure_axis_tick_size);
xlim([min(xi) max(xi)])
ylim([0 0.065])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tighten up the figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig1 = tightfig;
fig1.Position = FigureSize;
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];
fig1.PaperPositionMode = 'auto';

if shouldSaveFigures == 1
print('-depsc2', 'figures/tfit_error.eps')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Error distance distribution figures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FigureSize = [50 50 figure_width_1col+4 150*scaleFactor];

fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error PDF plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r=linspace(0,100,1000);
pdf1d = tdistpdf(r/sigma_t,nu)/sigma_t;
r(end+1)=1000; pdf1d(end+1) = 0.0;
t2d_pdf = @(z) interp1(r,pdf1d,z);

cdf_2d = cumtrapz(r,pdf1d);
t2d_cdf = @(z) interp1(r,cdf_2d,z);
cdf_2d(end+1) = 1.0; % so the interpolation algorithm has something to hang its hat on.

data =  dist_error_despiked;
histwidth = 50;
nbins = 100;

% Create the bins for the data
binwidth = histwidth/nbins;
edges = [(0:binwidth:(histwidth-binwidth))';histwidth*10];
binleft = edges(1:end-1);

% plot the data
count = histcounts(data,edges)';
g = bar(binleft, count/(length(data)*binwidth), 'histc'); hold on;
g.FaceColor = 0.8*[1.0 1.0 1.0];

% create bins for the analytical pdf
xi_mid = linspace(0,histwidth-binwidth,100)';
xi_right = linspace(histwidth-binwidth,histwidth,10)';
xi = [xi_mid;xi_right];

% plot the analytical pdf
pdf = t2d_pdf;
edgedensity = integral(pdf,(histwidth-binwidth),20*histwidth)/binwidth;
denfunc = edgedensity*ones(size(xi));
denfunc(1:100) = pdf(xi_mid);
plot(xi,denfunc,'LineWidth',2*scaleFactor,'Color',0.0*[1.0 1.0 1.0])

v95 = interp1(t2d_cdf(r),r, 0.95);
plot([v95 v95],get(gca,'ylim'), 'LineWidth', 0.5*scaleFactor, 'Color', 0.4*[1.0 1.0 1.0]);

xlabel('meters', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
set( gca, 'FontSize', figure_axis_tick_size);

fig1 = tightfig;
fig1.Position = FigureSize;
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];
fig1.PaperPositionMode = 'auto';

if shouldSaveFigures == 1
print('-depsc2', 'figures/tfit_distance_error.eps')
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Position error histogram
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

subplot(2,2,1)
plot_hist_with_pdf( error, t_pdf, 60, 100 )
% plot_hist_with_pdf( error_big, gaussian_pdf_big, 60, 100 )

subplot(2,2,2)
plot_hist_with_pdf( error_small, position_pdf_small, 60, 100 )
% plot_hist_with_pdf( error_small, gaussian_pdf_small, 20, 100 )

subplot(2,2,3)
plot_hist_with_pdf( a_big, velocity_pdf, 10e-5, 50 )

subplot(2,2,4)
plot_hist_with_pdf( a_small, velocity_pdf_small, 10e-5, 50 )
% plot_hist_with_pdf( a_small, exponential_pdf_small, 10e-5, 50 )


figure
subplot(1,2,1)
plot_hist_with_cdf( a_big, velocity_cdf, 10e-5, 50 )

subplot(1,2,2)
plot_hist_with_cdf( a_small, velocity_cdf_small, 10e-5, 50 )
