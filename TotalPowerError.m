function totalError = TotalPowerError( sigma, nu, a, T, S, drifters, shouldDisplay)

a = 10^(a);
w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));
Ndrifters = length(drifters.x);

drifters.uSmoothedTotalVariance = zeros(Ndrifters,1);
drifters.vSmoothedTotalVariance = zeros(Ndrifters,1);
drifters.axSmoothedTotalVariance = zeros(Ndrifters,1);
drifters.aySmoothedTotalVariance = zeros(Ndrifters,1);
for iDrifter = 1:Ndrifters
    x = drifters.x{iDrifter};
    y = drifters.y{iDrifter};
    t = drifters.t{iDrifter};
    
    tension = zeros(S,1);
    tension(T) = 1/a^2;
    [m_x,m_y,~,~,B,~,~] = bspline_bivariate_fit_with_tension(t,x,y,ones(size(x))*sigma,ones(size(x))*sigma,S,tension, w);
    
    X = squeeze(B(:,:,1));
    x2 = X*m_x; y2 = X*m_y;
    
    D = FiniteDifferenceMatrixNoBoundary(1,t,1);
    drifters.uSmoothedTotalVariance(iDrifter) = sum((D*x2).^2);
    drifters.vSmoothedTotalVariance(iDrifter) = sum((D*y2).^2);
    
    D = FiniteDifferenceMatrixNoBoundary(2,t,1);
    drifters.axSmoothedTotalVariance(iDrifter) = sum((D*x2).^2);
    drifters.aySmoothedTotalVariance(iDrifter) = sum((D*y2).^2);
end

ratio_u = (drifters.uTotalVariance - drifters.uSmoothedTotalVariance)./drifters.expectedTotalVelocityNoiseVariance;
ratio_v = (drifters.vTotalVariance - drifters.vSmoothedTotalVariance)./drifters.expectedTotalVelocityNoiseVariance;
ratio_ax = (drifters.axTotalVariance - drifters.axSmoothedTotalVariance)./drifters.expectedTotalAccelerationNoiseVariance;
ratio_ay = (drifters.ayTotalVariance - drifters.aySmoothedTotalVariance)./drifters.expectedTotalAccelerationNoiseVariance;

lambda = mean([ratio_u; ratio_v; ratio_ax; ratio_ay]);

if shouldDisplay == 1
    fprintf('\t(acceleration, lambda) = (%g, %f)\n', a, lambda);
end


totalError = abs(log10(lambda));

