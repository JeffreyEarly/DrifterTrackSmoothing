% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults
addpath('support')

% Drifter to highlight in the final plots
choiceDrifter = 6;

shouldSaveFigures = 1;

drifters = load('sample_data/rho1_drifters_projected_ungridded.mat');

S = 3; % order of the spline
K = S+1;
T = 2; % order of the tension
nu = 5.5; sigma =  8;
variance_of_the_noise = sigma*sigma*nu/(nu-2);
w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));

% Pull out the data of interest
x = drifters.x{choiceDrifter};
y = drifters.y{choiceDrifter};
t = drifters.t{choiceDrifter};
N = length(t);

% Estimate the velocities...
[u_rms_estimate_spectral, u_std, u_mean] = EstimateRMSVelocityFromSpectrum(t,x,sqrt(variance_of_the_noise), 0);
[v_rms_estimate_spectral, v_std, v_mean] = EstimateRMSVelocityFromSpectrum(t,y,sqrt(variance_of_the_noise), 0);

% ...and accelerations
[ax_rms_estimate_spectral, ax_std, ax_mean] = EstimateRMSAccelerationFromSpectrum(t,x,sqrt(variance_of_the_noise));
[ay_rms_estimate_spectral, ay_std, ay_mean] = EstimateRMSAccelerationFromSpectrum(t,y,sqrt(variance_of_the_noise));

sigma3Threshold = 1-0.9973;
gps_cdf = @(z) abs(tcdf(z/sigma,nu) - sigma3Threshold/2);
sigma3equiv = 3*sigma; % this would be correct if it were Gaussian
sigma3equiv = abs(fminsearch( gps_cdf, -50, optimset('TolX', 0.001, 'TolFun', 0.001) ));

% % Now estimate the optimal tension parameter
% % NOTE -- this is for outliers, so we use a large expected DOF.
% dt = median(diff(t));
% expectedDOFx = 1 + sigma3equiv/(u_rms_estimate_spectral*dt);
% lambda_x = 1/(ax_rms_estimate_spectral^2);
% 
% expectedDOFy = 1 + sigma3equiv/(v_rms_estimate_spectral*dt);
% lambda_y = 1/(ay_rms_estimate_spectral^2);
% 
% % Now we set a threshold for what constitutes an outlier. In this case we
% % choose points that have 1 in 10000 odds of occurring.
outlierThreshold = 0.0001;
gps_cdf = @(z) abs(tcdf(z/sigma,nu) - outlierThreshold/2);
range(1) = fminsearch( gps_cdf, -50, optimset('TolX', 0.001, 'TolFun', 0.001) );
range(2) = -range(1); % it's symmetric
 
 lambda_x = 7.2947e+10; % value determined after optimization
% errorFunction = @(log10lambda) KolmogorovSmirnovErrorForTDistribution(sigma,nu,log10lambda,T,S,range,drifters.t,drifters.x,1);
% optimalLog10lambda = fminsearch( errorFunction, log10(lambda_x), optimset('TolX', 0.01, 'TolFun', 0.01) );
% lambda_x = 10^(optimalLog10lambda);
% 
 lambda_y = 5.4806e+10; % value determined after optimization
% errorFunction = @(log10lambda) KolmogorovSmirnovErrorForTDistribution(sigma,nu,log10lambda,T,S,range,drifters.t,drifters.y,1);
% optimalLog10lambda = fminsearch( errorFunction, log10(lambda_y), optimset('TolX', 0.01, 'TolFun', 0.01) );
% lambda_y = 10^(optimalLog10lambda);

epsilon_high_tension = [];
epsilon = [];
Ndrifters = length(drifters.x);
x_interp = cell(Ndrifters,1);
y_interp = cell(Ndrifters,1);
u_interp = cell(Ndrifters,1);
v_interp = cell(Ndrifters,1);
ax_interp = cell(Ndrifters,1);
ay_interp = cell(Ndrifters,1);
t_interp = cell(Ndrifters,1);
rejectedPointIndices_x = cell(Ndrifters,1);
rejectedPointIndices_y = cell(Ndrifters,1);
for iDrifter = 1:Ndrifters
    x = drifters.x{iDrifter};
    y = drifters.y{iDrifter};
    t = drifters.t{iDrifter};
    N = length(t);

    t_knot = NaturalKnotsForSpline( t, K, 1 );
    [m_x,Cm_x,Bx,Bxq,txq] = bspline_fit_with_tension(t,x,sigma,t_knot,S,T,lambda_x,w);
    [m_y,Cm_y,By,Byq,tyq] = bspline_fit_with_tension(t,y,sigma,t_knot,S,T,lambda_y,w);

    X = squeeze(Bx(:,:,1));
    Y = squeeze(By(:,:,1));

    epsilon_x = x-X*m_x;
    epsilon_y = y-Y*m_y;
    epsilon_high_tension = [epsilon_high_tension; epsilon_x; epsilon_y];
    
    % Reject interior points only--leave endpoints untouched.
    rejectedPointIndices_x{iDrifter} = find(epsilon_x(2:end-1) < range(1) | epsilon_x(2:end-1) > range(2) );
    rejectedPointIndices_y{iDrifter} = find(epsilon_y(2:end-1) < range(1) | epsilon_y(2:end-1) > range(2) );
    if ~isempty(rejectedPointIndices_x{iDrifter})
        rejectedPointIndices_x{iDrifter} = rejectedPointIndices_x{iDrifter}+1;
    end
    if ~isempty(rejectedPointIndices_y{iDrifter})
        rejectedPointIndices_y{iDrifter} = rejectedPointIndices_y{iDrifter}+1;
    end
    
    t_x = t; t_x(rejectedPointIndices_x{iDrifter}) = [];
    x_reduced = x; x_reduced(rejectedPointIndices_x{iDrifter}) = [];
    t_y = t; t_y(rejectedPointIndices_y{iDrifter}) = [];
    y_reduced = y; y_reduced(rejectedPointIndices_y{iDrifter}) = [];
    
    % Estimate the velocities...
    u_estimate_spectral = EstimateRMSVelocityFromSpectrum(t_x,x_reduced,sqrt(variance_of_the_noise), 0);
    v_estimate_spectral = EstimateRMSVelocityFromSpectrum(t_y,y_reduced,sqrt(variance_of_the_noise), 0);

    % ...and accelerations
    [ax_rms_estimate_spectral, ax_std, ax_mean] = EstimateRMSAccelerationFromSpectrum(t_x,x_reduced,sqrt(variance_of_the_noise));
    [ay_rms_estimate_spectral, ay_std, ay_mean] = EstimateRMSAccelerationFromSpectrum(t_y,y_reduced,sqrt(variance_of_the_noise));
    
    dt_rms = sqrt(mean(diff(t_x).^2));
    expectedDOFx = 1 + sigma3equiv/(u_estimate_spectral*dt_rms);
    lambda_x2 = (expectedDOFx-1)/(expectedDOFx*ax_std^2);

    dt_rms = sqrt(mean(diff(t_y).^2));
    expectedDOFy = 1 + sigma3equiv/(u_estimate_spectral*dt_rms);
    lambda_y2 = (expectedDOFy-1)/(expectedDOFy*ay_std^2);
    
    % Now create a grid that will coincide for all drifters, using the fact
    % that zero coincides for all them. But also include the end points.
    res = 5*60;
    tq = res*( ceil(min(t_x(1),t_y(2))/res):1:floor(max(t_x(end),t_y(end))/res) )';
    if min(t) < min(tq)
        tq = [min(t); tq];
    end
    if max(t) > max(tq)
        tq(end+1) = max(t);
    end    
    
    t_interp{iDrifter} = tq;
    
    % Fit with the new reduced tension
    Sigma_x = variance_of_the_noise*eye(length(t_x));
    t_knot_x = NaturalKnotsForSpline( t_x, K, 1 );
    [m_x,Cm_x,Bx,~,~,Wx] = bspline_fit_with_tension(t_x,x_reduced,sigma,t_knot_x,S,T,lambda_x2,w,ax_mean);
    X = squeeze(Bx(:,:,1));
    SE_x = X*Cm_x*X.'*Wx*Sigma_x;
    dof_x = variance_of_the_noise/ mean(diag(SE_x));
    Bq = bspline(tq,t_knot_x,K);
    Xq = squeeze(Bq(:,:,1));
    x_interp{iDrifter} = Xq*m_x;
    u_interp{iDrifter} = squeeze(Bq(:,:,2))*m_x;
    ax_interp{iDrifter} = squeeze(Bq(:,:,3))*m_x;
    
    Sigma_y = variance_of_the_noise*eye(length(t_y));
    t_knot_y = NaturalKnotsForSpline( t_y, K, 1 );
    [m_y,Cm_y,By,~,~,Wy] = bspline_fit_with_tension(t_y,y_reduced,sigma,t_knot_y,S,T,lambda_y2,w,ay_mean);
    Y = squeeze(By(:,:,1));
    SE_y = Y*Cm_y*Y.'*Wy*Sigma_y;
    dof_y = variance_of_the_noise/mean(diag(SE_y));
    Bq = bspline(tq,t_knot_y,K);
    Yq = squeeze(Bq(:,:,1));
    y_interp{iDrifter} = Yq*m_y;
    v_interp{iDrifter} = squeeze(Bq(:,:,2))*m_y;
    ay_interp{iDrifter} = squeeze(Bq(:,:,3))*m_y;
    
    
    epsilon_x = x_reduced-X*m_x;
    epsilon_y = y_reduced-Y*m_y;
    epsilon = [epsilon; epsilon_x; epsilon_y];
end

ax = [];
ay = [];
for iDrifter = 1:Ndrifters
   ax = [ax;ax_interp{iDrifter}];
   ay = [ay;ay_interp{iDrifter}];
end

t_pdf = @(z) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2));

s = 1/1000; % scale
for iDrifter = 1:Ndrifters
figure
plot(s*x_interp{iDrifter},s*y_interp{iDrifter}, 'LineWidth', 0.5*scaleFactor, 'Color',0.4*[1.0 1.0 1.0]), hold on
scatter(s*drifters.x{iDrifter},s*drifters.y{iDrifter},(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
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
plot(t_interp{choiceDrifter}/3600,s*x_interp{choiceDrifter}, 'LineWidth', 0.5*scaleFactor, 'Color',0.4*[1.0 1.0 1.0]), hold on
scatter(drifters.t{choiceDrifter}(rejectedPointIndices_x{choiceDrifter})/3600,s*drifters.x{choiceDrifter}(rejectedPointIndices_x{choiceDrifter}),(6.5*scaleFactor)^2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
scatter(drifters.t{choiceDrifter}/3600,s*drifters.x{choiceDrifter},(2.5*scaleFactor)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
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

data = epsilon;
histwidth = 80;
nbins = 100;

% Create the bins for the data
binwidth = histwidth/nbins;
edges = [-histwidth*100;((-histwidth/2+binwidth):binwidth:(histwidth/2-binwidth))';histwidth*100];
binleft = linspace((-histwidth/2),(histwidth/2-binwidth),nbins)';

% plot the data

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
ylim([0 0.175])