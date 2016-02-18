addpath('./support');

drifters = load('sample_data/rho1_drifters_projected_ungridded.mat');
output = sprintf('rho1_drifters_smoothed_interpolated.mat');

Ndrifters = length(drifters.x);

% Tension order
T = 2; a_start = log10(4.26e-6);
% T = 3; a_start = log10(2.99e-9);
% T = 4; a_start = log10(1.96e-12);
S = T+1; K = S+1;
nu = 5.5; sigma =  8.0;

% a = 10^a_start;
% outlierCut = 93;

w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));

% Find the point at which we'll reject 1 in 1000 points
fprintf('Generating the 2D student t-distribution...\n')
[r, pdf1d] = TwoDimStudentTProbabilityDistributionFunction( sigma, nu, 150, 3001 );
cdf_2d = cumtrapz(r,pdf1d);
r(end+1)=20000; cdf_2d(end+1) = 1.0; % so the interpolation algorithm has something to hang its hat on.
outlierCut = interp1(cdf_2d,r, 0.9999);

% Optimize the tension parameter so we can get rid of the outliers
fprintf('Seaching for the optimal tension parameter using distances less than %f...\n', outlierCut)
errorFunction = @(a) KolmogorovSmirnovErrorFor2DTDistribution( sigma, nu, a, T, S, outlierCut, drifters, r, cdf_2d, 1);
optimalAcceleration = fminsearch( errorFunction, a_start, optimset('TolX', 0.1, 'TolFun', 0.1) );
fprintf('Optimal acceleration tension is %g\n', 10^(optimalAcceleration(1)) );
a = 10^(optimalAcceleration(1));

fprintf('Removing points outside of our cutoff...\n')
despikedDrifters.x = cell(Ndrifters,1);
despikedDrifters.y = cell(Ndrifters,1);
despikedDrifters.t = cell(Ndrifters,1);
despikedDrifters.xTotalVariance = zeros(Ndrifters,1);
despikedDrifters.yTotalVariance = zeros(Ndrifters,1);
despikedDrifters.expectedTotalPositionNoiseVariance = zeros(Ndrifters,1);

despikedDrifters.uTotalVariance = zeros(Ndrifters,1);
despikedDrifters.vTotalVariance = zeros(Ndrifters,1);
despikedDrifters.expectedTotalVelocityNoiseVariance = zeros(Ndrifters,1);

despikedDrifters.axTotalVariance = zeros(Ndrifters,1);
despikedDrifters.ayTotalVariance = zeros(Ndrifters,1);
despikedDrifters.expectedTotalAccelerationNoiseVariance = zeros(Ndrifters,1);
for iDrifter = 1:Ndrifters
    x = drifters.x{iDrifter};
    y = drifters.y{iDrifter};
    t = drifters.t{iDrifter};
    
    tension = zeros(S,1);
    tension(T) = 1/a^2;
    [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = bspline_bivariate_fit_with_tension(t,x,y,ones(size(x))*sigma,ones(size(x))*sigma,S,tension, w);
    
    X = squeeze(B(:,:,1));
    error_x_big = X*m_x - x;
    error_y_big = X*m_y - y;
    dist_error = sqrt( error_x_big.*error_x_big + error_y_big.*error_y_big );
    goodDrifters = dist_error < outlierCut;
    
    despikedDrifters.x{iDrifter} = x(goodDrifters);
    despikedDrifters.y{iDrifter} = y(goodDrifters);
    despikedDrifters.t{iDrifter} = t(goodDrifters);
    
    x = despikedDrifters.x{iDrifter};
    y = despikedDrifters.y{iDrifter};
    t = despikedDrifters.t{iDrifter};
    
    dx = (sigma*sigma*nu/(nu-2))*ones(size(t));
    
    despikedDrifters.xTotalVariance(iDrifter) = sum(x.*x);
    despikedDrifters.yTotalVariance(iDrifter) = sum(y.*y);
    despikedDrifters.expectedTotalPositionNoiseVariance(iDrifter) = sum(dx);
    
    D = FiniteDifferenceMatrixNoBoundary(1,t,1);
    despikedDrifters.uTotalVariance(iDrifter) = sum((D*x).^2);
    despikedDrifters.vTotalVariance(iDrifter) = sum((D*y).^2);
    despikedDrifters.expectedTotalVelocityNoiseVariance(iDrifter) = sum((D.^2)*dx);
    
    D = FiniteDifferenceMatrixNoBoundary(2,t,1);
    despikedDrifters.axTotalVariance(iDrifter) = sum((D*x).^2);
    despikedDrifters.ayTotalVariance(iDrifter) = sum((D*y).^2);
    despikedDrifters.expectedTotalAccelerationNoiseVariance(iDrifter) = sum((D.^2)*dx);
    
    fprintf('Decreased points from %d to %d\n',length(drifters.t{iDrifter}),sum(goodDrifters))
end


fprintf('Seaching for a new optimal tension parameter to preserve total variance...\n')
errorFunction = @(a) TotalPowerError( sigma, nu, a, T, S, despikedDrifters, 1);
optimalAcceleration = fminsearch( errorFunction, log10(a), optimset('TolX', 0.1, 'TolFun', 0.1) );
fprintf('Optimal acceleration tension is %g\n', 10^(optimalAcceleration(1)) );
a = 10^(optimalAcceleration(1));


fprintf('Performing final fit and gridding the output...\n')
x_interp = cell(Ndrifters,1);
y_interp = cell(Ndrifters,1);
t_interp = cell(Ndrifters,1);
u = cell(Ndrifters,1);
v = cell(Ndrifters,1);
ax = cell(Ndrifters,1);
ay = cell(Ndrifters,1);
x_error = cell(Ndrifters,1);
y_error = cell(Ndrifters,1);
x_error_despiked = cell(Ndrifters,1);
y_error_despiked = cell(Ndrifters,1);

despikedDrifters.xSmoothedTotalVariance = zeros(Ndrifters,1);
despikedDrifters.ySmoothedTotalVariance = zeros(Ndrifters,1);

despikedDrifters.uSmoothedTotalVariance = zeros(Ndrifters,1);
despikedDrifters.vSmoothedTotalVariance = zeros(Ndrifters,1);

despikedDrifters.axSmoothedTotalVariance = zeros(Ndrifters,1);
despikedDrifters.aySmoothedTotalVariance = zeros(Ndrifters,1);

for iDrifter = 1:Ndrifters
    x = despikedDrifters.x{iDrifter};
    y = despikedDrifters.y{iDrifter};
    t = despikedDrifters.t{iDrifter};
    N = length(t);
    
    dx = ones(size(x))*sigma;
    dy = ones(size(y))*sigma;
    tension = zeros(S,1);
    tension(T) = 1/a^2;
    [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = bspline_bivariate_fit_with_tension(t,x,y,dx,dy,S,0.1*tension, w);
    
    X = squeeze(B(:,:,1));
    x_error_despiked{iDrifter} = X*m_x - x;
    y_error_despiked{iDrifter} = X*m_y - y;
    
    x2 = X*m_x; y2 = X*m_y;
    despikedDrifters.xSmoothedTotalVariance(iDrifter) = sum(x2.*x2);
    despikedDrifters.ySmoothedTotalVariance(iDrifter) = sum(y2.*y2);
    
    D = FiniteDifferenceMatrixNoBoundary(1,t,1);
    despikedDrifters.uSmoothedTotalVariance(iDrifter) = sum((D*x2).^2);
    despikedDrifters.vSmoothedTotalVariance(iDrifter) = sum((D*y2).^2);
    
    D = FiniteDifferenceMatrixNoBoundary(2,t,1);
    despikedDrifters.axSmoothedTotalVariance(iDrifter) = sum((D*x2).^2);
    despikedDrifters.aySmoothedTotalVariance(iDrifter) = sum((D*y2).^2);
    
    % Recreate that knots that were used internally
    t_knot = NaturalKnotsForSpline( t, K );
    
    % Now locate at ALL the errors, not just the de-spiked ones
    B = bspline(drifters.t{iDrifter},t_knot,K);
    X = squeeze(B(:,:,1));
    x_error{iDrifter} = X*m_x - drifters.x{iDrifter};
    y_error{iDrifter} = X*m_y - drifters.y{iDrifter};
    
    % Now create a grid that will coincide for all drifters, using the fact
    % that zero coincides for all them. But also include the end points.
    res = 5*60;
    tq = res*( ceil(min(t)/res):1:floor(max(t)/res) )';
    if min(t) < min(tq)
        tq = [min(t); tq];
    end
    if max(t) > max(tq)
        tq(end+1) = max(t);
    end
    Bq = bspline(tq,t_knot,K);
    
    Xq = squeeze(Bq(:,:,1));
    x_interp{iDrifter} = Xq*m_x;
    y_interp{iDrifter} = Xq*m_y;
    t_interp{iDrifter} = tq;
    
    Vq = squeeze(Bq(:,:,2));
    u{iDrifter} = Vq*m_x;
    v{iDrifter} = Vq*m_y;
    
    Aq = squeeze(Bq(:,:,3));
    ax{iDrifter} = Aq*m_x;
    ay{iDrifter} = Aq*m_y;
end

ratio_x = (despikedDrifters.xTotalVariance - despikedDrifters.xSmoothedTotalVariance)./despikedDrifters.expectedTotalPositionNoiseVariance;
ratio_u = (despikedDrifters.uTotalVariance - despikedDrifters.uSmoothedTotalVariance)./despikedDrifters.expectedTotalVelocityNoiseVariance;
ratio_v = (despikedDrifters.vTotalVariance - despikedDrifters.vSmoothedTotalVariance)./despikedDrifters.expectedTotalVelocityNoiseVariance;
ratio_ax = (despikedDrifters.axTotalVariance - despikedDrifters.axSmoothedTotalVariance)./despikedDrifters.expectedTotalAccelerationNoiseVariance;
ratio_ay = (despikedDrifters.ayTotalVariance - despikedDrifters.aySmoothedTotalVariance)./despikedDrifters.expectedTotalAccelerationNoiseVariance;

x_raw = drifters.x;
y_raw = drifters.y;
t_raw = drifters.t;

x = x_interp;
y = y_interp;
t = t_interp;

f0 = drifters.f0;
lastDeployment = drifters.lastDeployment;
lat0 = drifters.lat0;
lon0 = drifters.lon0;
maxExperimentLength = drifters.maxExperimentLength;

save(output,'f0','lat0','lon0','maxExperimentLength', 't', 'x', 'y', 'u', 'v', 'ax', 'ay', 'x_raw', 'y_raw', 't_raw', 'x_error', 'y_error', 'x_error_despiked', 'y_error_despiked', 'S', 'nu', 'sigma', 'T', 'a', 'outlierCut')