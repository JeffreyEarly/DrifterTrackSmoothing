function totalError = GaussianFitScalarFunction( sigma, a, t,x,y,S,errorMetric)

w = @(z)(sigma*sigma);
gamma = [0; 1/a^2; 0];    
dx = ones(size(x))*sigma;
dy = ones(size(y))*sigma;

[m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = bspline_bivariate_fit_with_tension(t,x,y,dx,dy,S,gamma, w);

a_big = [squeeze(Bq(:,:,3))*m_x; squeeze(Bq(:,:,3))*m_y];

X = squeeze(B(:,:,1));
error_x_big = X*m_x - x;
error_y_big = X*m_y - y;
error_big = [error_x_big;error_y_big];

if (strcmp(errorMetric,'total'))
    totalError = abs(log10(a/std(a_big))) + abs(log10(sigma/std(error_big)));
elseif (strcmp(errorMetric,'acceleration'))
    totalError = abs(log10(a/std(a_big)));
else
    totalError = abs(log10(sigma/std(error_big)));
end