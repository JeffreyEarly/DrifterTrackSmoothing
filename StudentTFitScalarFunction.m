function totalError = StudentTFitScalarFunction( sigma, nu, a, drifters)

a = 10^(a);
S = 3; % order of the spline
maxlag = 15;
w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));
Ndrifters = length(drifters.x);

ACx_big = zeros(maxlag+1,1);
ACy_big = zeros(maxlag+1,1);
n = 0;
for iDrifter = 1:Ndrifters
    x = drifters.x{iDrifter};
    y = drifters.y{iDrifter};
    t = drifters.t{iDrifter};
    
    dx = ones(size(x))*sigma;
    dy = ones(size(y))*sigma;
    [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = bspline_bivariate_fit_with_tension(t,x,y,dx,dy,S,[0; 1/a^2; 0], w);
    
     X = squeeze(B(:,:,1));
    error_x_big = X*m_x - x;
    error_y_big = X*m_y - y;
    
    ACx_big = ACx_big + Autocorrelation(error_x_big,maxlag);
    ACy_big = ACy_big + Autocorrelation(error_y_big,maxlag);
    
    n = n + length(drifters.x{iDrifter});
end

ACx_big = ACx_big/Ndrifters;
ACy_big = ACy_big/Ndrifters;
AC = (ACx_big + ACy_big)/2;

[p, Q] = LjungBoxTest( AC, n );

totalError = Q(10);