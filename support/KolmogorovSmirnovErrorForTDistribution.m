% Return a measure of the Kolmogorov-Smirnov test metric for the
% distribution matched to the position errors of drifters, given a
% particular tension, a.
%
% the range allows you to restrict the range (in meters) over which the
% test is applied.
function totalError = KolmogorovSmirnovErrorForTDistribution( sigma, nu, log10lambda, T, S, range, t, x, shouldDisplay)

lambda = 10^(log10lambda);
w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));
gps_cdf = @(z) tcdf(z/sigma,nu);

Ndrifters = length(x);
error = [];
for iDrifter = 1:Ndrifters
    t_knot = NaturalKnotsForSpline( t{iDrifter}, S+1, 1 );
    [m_x,~,B,~,~] = bspline_fit_with_tension(t{iDrifter},x{iDrifter},sigma,t_knot,S,T,lambda,w);
    error = [error; squeeze(B(:,:,1))*m_x - x{iDrifter}];
end

x = sort(error( error > range(1) & error < range(2) ));
n = length(x);
y_data = (1:n)'/n;
y = gps_cdf(x);

D = max(abs(y-y_data));

totalError = (sqrt(n) + 0.12 + 0.11/sqrt(n))*D;

if shouldDisplay == 1
    fprintf('(acceleration, totalError) = (%g, %f)\n', log10lambda, totalError);
end
