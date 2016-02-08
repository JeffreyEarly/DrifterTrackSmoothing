% Using the Ljung-Box test statistic

% nu = 2.00; sigma = 8;
% errorFunction = @(a) StudentTFitScalarFunction( sigma, nu, a, drifters);
% optimalAcceleration = fminsearch( errorFunction, log10(0.9e-5) )
% a = 1.1994e-05;

% nu = 2.00; sigma = 4;
% errorFunction = @(a) StudentTFitScalarFunction( sigma, nu, a, drifters);
% optimalAcceleration = fminsearch( errorFunction, log10(1.7e-5) )
% a= 10^(-4.9666)

% nu = 2.00; sigma = 1.7;
% errorFunction = @(a) StudentTFitScalarFunction( sigma, nu, a, drifters);
% optimalAcceleration = fminsearch( errorFunction, log10(1.0e-5) )
% a = 9.3254e-06;

nu = 2.00; sigma = 1.25; a = 7.2975e-06;
errorFunction = @(a) StudentTFitScalarFunction( sigma, nu, a, drifters);
optimalAcceleration = fminsearch( errorFunction, log10(0.9e-5) )

% nu = 2.00; sigma = 2.5; a = 8.9794e-06;
% errorFunction = @(a) StudentTFitScalarFunction( sigma, nu, a, drifters);
% optimalAcceleration = fminsearch( errorFunction, log10(1.0e-5) )


% nu = 2.00; sigma = 5.0; a = 1.1272e-05;
% errorFunction = @(a) StudentTFitScalarFunction( sigma, nu, a, drifters);
% optimalAcceleration = fminsearch( errorFunction, log10(8e-6) )
% 
% 
% nu = 2.00; sigma = 10.0; a = 1.2317e-05;
% errorFunction = @(a) StudentTFitScalarFunction( sigma, nu, a, drifters);
% optimalAcceleration = fminsearch( errorFunction, log10(7e-6) )
