shouldReoptimizeAfterDespiking = 0;

drifters = load('sample_data/projected_ungridded_rho1_drifters.mat');

Ndrifters = length(drifters.x);

% Tension order
T = 2; a_start = log10(4.26e-6);
S = T; K = S+1;
nu = 5.5; sigma =  1.0;

targetAcceleration = 3e-5;
outlierCut = 300;

w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));

% fprintf('Seaching for the optimal tension parameter using target acceleration %f...\n', targetAcceleration)
% errorFunction = @(a) AccelerationErrorForDespiking( sigma, nu, a, T, S, targetAcceleration, drifters, 1);
% optimalAcceleration = fminsearch( errorFunction, a_start, optimset('TolX', 0.1, 'TolFun', 0.1) );
% fprintf('Optimal acceleration tension is %g\n', 10^(optimalAcceleration(1)) );
% a = 10^(optimalAcceleration(1));

a = 3.12712e-06;

fprintf('Removing points outside of our cutoff...\n')
despikedDrifters.x = cell(Ndrifters,1);
despikedDrifters.y = cell(Ndrifters,1);
despikedDrifters.t = cell(Ndrifters,1);
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
    
    fprintf('Decreased points from %d to %d\n',length(t),sum(goodDrifters))
end

if shouldReoptimizeAfterDespiking == 1
    fprintf('Seaching for a new optimal tension parameter...\n')
    errorFunction = @(a) AccelerationErrorForDespiking( sigma, nu, a, T, S, targetAcceleration, drifters, 1);
    optimalAcceleration = fminsearch( errorFunction, log10(a), optimset('TolX', 0.1, 'TolFun', 0.1) );
    fprintf('Optimal acceleration tension is %g\n', 10^(optimalAcceleration(1)) );
    a = 10^(optimalAcceleration(1));
end

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
for iDrifter = 1:Ndrifters
    x = despikedDrifters.x{iDrifter};
    y = despikedDrifters.y{iDrifter};
    t = despikedDrifters.t{iDrifter};
    N = length(t);
    
    dx = ones(size(x))*sigma;
    dy = ones(size(y))*sigma;
    tension = zeros(S,1);
    tension(T) = 0.01*1/a^2;
    [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = bspline_bivariate_fit_with_tension(t,x,y,dx,dy,S,tension, w);
    
    X = squeeze(B(:,:,1));
    x_error_despiked{iDrifter} = X*m_x - x;
    y_error_despiked{iDrifter} = X*m_y - y;
    
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

save(sprintf('smoothed_interpolated_rho1_drifters_DespikedLowTension.mat'),'f0','lat0','lon0','maxExperimentLength', 't', 'x', 'y', 'u', 'v', 'ax', 'ay', 'x_raw', 'y_raw', 't_raw', 'x_error', 'y_error', 'x_error_despiked', 'y_error_despiked', 'S', 'nu', 'sigma', 'T', 'a', 'outlierCut')