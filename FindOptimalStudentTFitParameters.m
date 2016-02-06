% Using the Ljung-Box test statistic

% nu = 2.00; sigma = 8;
% errorFunction = @(a) StudentTFitScalarFunction( sigma, nu, a, drifters);
% optimalAcceleration = fminsearch( errorFunction, log10(0.9e-5) )
% a = 1.1994e-05;

% nu = 2.00; sigma = 4;
% errorFunction = @(a) StudentTFitScalarFunction( sigma, nu, a, drifters);
% optimalAcceleration = fminsearch( errorFunction, log10(1.7e-5) )
% a= 10^(-4.9666)

nu = 2.00; sigma = 1.7;
errorFunction = @(a) StudentTFitScalarFunction( sigma, nu, a, drifters);
optimalAcceleration = fminsearch( errorFunction, log10(1.0e-5) )
% a = 9.3254e-06;