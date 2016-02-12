% Return a measure of the Kolmogorov-Smirnov test metric for the
% distribution matched to the position errors of drifters, given a
% particular tension, a.
%
% the range allows you to restrict the range (in meters) over which the
% test is applied.
function totalError = KolmogorovSmirnovErrorFor2DTDistribution( sigma, nu, a, cut, drifters, r, cdf, shouldDisplay)

a = 10^(a);
S = 3; % order of the spline
w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));
Ndrifters = length(drifters.x);

dist_error = [];
for iDrifter = 1:Ndrifters
    x = drifters.x{iDrifter};
    y = drifters.y{iDrifter};
    t = drifters.t{iDrifter};
    
    tension = zeros(S,1);
    tension(2) = 1/a^2;
    [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = bspline_bivariate_fit_with_tension(t,x,y,ones(size(x))*sigma,ones(size(x))*sigma,S,tension, w);
    
    X = squeeze(B(:,:,1));
    error_x_big = X*m_x - x;
    error_y_big = X*m_y - y;
    dist_error = [dist_error; sqrt( error_x_big.*error_x_big + error_y_big.*error_y_big )];
end

x = sort(dist_error( dist_error < cut ));
n = length(x);
y_data = (1:n)'/n;
y = interp1(r,cdf,x);
D = max(abs(y-y_data));

lambda = (sqrt(n) + 0.12 + 0.11/sqrt(n))*D;

if shouldDisplay == 1
    fprintf('\t(acceleration, lambda) = (%g, %f)\n', a, lambda);
end

% j=1:25;
% p = sum(2*((-1).^(j-1)).*exp(-2*j.*j*lambda*lambda))

totalError = lambda;