% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults
addpath('support')

% Drifter to highlight in the final plots
choiceDrifter = 6;

shouldSaveFigures = 0;

% How many data points do we have total
drifters = load('sample_data/projected_ungridded_rho1_drifters.mat');

T = 2; a_start = log10(0.4e-5);
% T = 3; a_start = log10(2.99e-9);
T = 4; a_start = log10(1.96e-12);
S = T+1; K = S+1;
nu = 5.5; sigma = 8.0;

w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));
a = 10^a_start;
a = 4e-12;

NDrifters = 9;
maxT = [];
for iDrifter = 1:NDrifters
    maxT(end+1) = max(drifters.t{iDrifter});
end
lastTime = min(maxT);


% Find the point at which we'll reject 1 in 1000 points
fprintf('Generating the 2D student t-distribution...\n')
[r, pdf1d] = TwoDimStudentTProbabilityDistributionFunction( sigma, nu, 150, 3001 );
cdf_2d = cumtrapz(r,pdf1d);
r(end+1)=20000; cdf_2d(end+1) = 1.0; % so the interpolation algorithm has something to hang its hat on.
outlierCut = interp1(cdf_2d,r, 0.9999);

% fprintf('Seaching for the optimal tension parameter using distances less than %f...\n', outlierCut)
% errorFunction = @(a) KolmogorovSmirnovErrorFor2DTDistribution( sigma, nu, a, T, S, outlierCut, drifters, r, cdf_2d, 1);
% optimalAcceleration = fminsearch( errorFunction, a_start, optimset('TolFun', 1) );
% fprintf('Optimal acceleration tension is %g\n', 10^(optimalAcceleration(1)) );
% a = 10^(optimalAcceleration(1));

% a = 4e-5;

res = 30*60;
tq = (0:res:lastTime)';
x = zeros(length(tq),NDrifters);
y = zeros(length(tq),NDrifters);
for iDrifter = 1:NDrifters

    d = ones(size(drifters.x{iDrifter}))*sigma;
    tension = zeros(S,1);
    tension(T) = 1/a^2;
    [m_x,m_y,Cm_x,Cm_y,B,~,~] = bspline_bivariate_fit_with_tension(drifters.t{iDrifter},drifters.x{iDrifter},drifters.y{iDrifter},d,d,S,tension, w);
    
    % Recreate that knots that were used internally
    t_knot = NaturalKnotsForSpline( drifters.t{iDrifter}, K );
    
    % Now create a grid that will coincide for all drifters, using the fact
    % that zero coincides for all them. But also include the end points.
    Bq = bspline(tq,t_knot,K);
    
    x(:,iDrifter) = squeeze(Bq(:,:,1))*m_x;
    y(:,iDrifter) = squeeze(Bq(:,:,1))*m_y;
end

[x_com, y_com, q, r] = CenterOfMass( x, y );

dt = tq(2)-tq(1);
cv = (diff(q,1,1) + sqrt(-1)*diff(r,1,1))/dt;

taper_bandwidth = 3;
[psi,lambda]=sleptap(size(cv,1),taper_bandwidth);
[omega_p,spp,snn,spn]=mspec(dt,cv,psi);

% convert from radians/second to cycles/day
f_p=omega_p*86400/(2*pi);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Spectrum
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure('Units', 'points', 'Position', [50 50 figure_width_2col 300*scaleFactor])
    set(gcf,'PaperPositionMode','auto')
    set(gcf, 'Color', 'w');
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Double the zero frequency for plotting purposes
	snn(1,:)=snn(1,:);
	spp(1,:)=spp(1,:);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    plot(f_p,spp, 'Color', [0.4 0.4 0.4], 'LineWidth', 2),ylog
	hold on
	plot(f_p,snn, 'Color', [0.4 0.4 0.4], 'LineWidth', 2)
	plot(f_p,vmean(spp,2), 'blue', 'LineWidth', 4)
	plot(f_p,vmean(snn,2), 'blue', 'LineWidth', 4)
    xlog
    
    xlim([min(f_p) max(f_p)])
    
    return

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
