drifters = load('sample_data/projected_ungridded_rho1_drifters.mat');

nu = 6; sigma_gps = 6; sigma = sigma_gps; S = 3;
errorFunction = @(a) StudentTFitMaximizedGPSFunction( sigma, nu, a, drifters);
optimalAcceleration = fminsearch( errorFunction, log10(0.4e-5), optimset('TolFun',1e-1) )
a = 5.3462e-06;

% Now optimizing against the T-distribution, but only at the core.
nu = 5; sigma_gps = 1.75; sigma = sigma_gps; S = 3;
errorFunction = @(a) StudentTFitMaximizedGPSFunction( sigma, nu, a, drifters);
optimalAcceleration = fminsearch( errorFunction, log10(0.4e-5), optimset('TolFun',1e-1) )
a = 3.9737e-06;

% Now optimizing against the T-distribution across all errors
nu = 5; sigma_gps = 1.75; sigma = sigma_gps; S = 3;
errorFunction = @(a) StudentTFitMaximizedGPSFunction( sigma, nu, a, drifters);
optimalAcceleration = fminsearch( errorFunction, log10(0.4e-5), optimset('TolFun',1e-1) )
a = 5.0192e-06

% Now optimizing against the T-distribution, but only at the core.
nu = 6; sigma_gps = 1.0; sigma = sigma_gps; S = 3;
errorFunction = @(a) StudentTFitMaximizedGPSFunction( sigma, nu, a, drifters);
optimalAcceleration = fminsearch( errorFunction, log10(0.4e-5), optimset('TolFun',1e-1) )
a = 2.2412e-06;

% drifters = load('sample_data/projected_ungridded_rho2_drifters.mat');
% nu = 10; sigma_gps = 1.75; sigma = sigma_gps; S = 3;
% errorFunction = @(a) StudentTFitMaximizedGPSFunction( sigma, nu, a, drifters);
% optimalAcceleration = fminsearch( errorFunction, log10(0.7e-5), optimset('TolFun',1e-2) )

% Using the Ljung-Box test statistic

% nu = 2.00; sigma = 8;
% errorFunction = @(a) StudentTFitScalarFunction( sigma, nu, a, drifters);
% optimalAcceleration = fminsearch( errorFunction, log10(0.9e-5) )
% a = 1.1994e-05;

% nu = 2.00; sigma = 4;
% errorFunction = @(a) StudentTFitScalarFunction( sigma, nu, a, drifters);
% % optimalAcceleration = fminsearch( errorFunction, log10(1.7e-5) )
% optimalAcceleration = fminsearch( errorFunction, log10(0.7e-5) )
% a= 10^(-4.9666)

% nu = 2.00; sigma = 4; sigma_gps = 1.75
% errorFunction = @(a) StudentTFitMaximizedGPSFunction( sigma, nu, a, drifters);
% optimalAcceleration = fminsearch( errorFunction, log10(0.7e-5), optimset('TolX',1e-2) )
% a = 9.2956e-06;

nu = 2.00; sigma = 4; sigma_gps = 4;
errorFunction = @(a) StudentTFitMaximizedGPSFunction( sigma, nu, a, drifters);
optimalAcceleration = fminsearch( errorFunction, log10(0.7e-5), optimset('TolX',1e-2) )
% a = 5.1614e-06

nu = 2.00; sigma_gps = 1.75; sigma = 3*sigma_gps;
errorFunction = @(a) StudentTFitMaximizedGPSFunction( sigma, nu, a, drifters);
optimalAcceleration = fminsearch( errorFunction, log10(0.7e-5), optimset('TolFun',1e-2) )
% a = 1.2913e-05;

nu = 5.00; sigma_gps = 1.75; sigma = sigma_gps;
errorFunction = @(a) StudentTFitMaximizedGPSFunction( sigma, nu, a, drifters);
optimalAcceleration = fminsearch( errorFunction, log10(0.7e-5), optimset('TolFun',1e-2) )
% a = 4.1839e-06;

nu = 3.5; sigma_gps = 1.75; sigma = sigma_gps;
errorFunction = @(a) StudentTFitMaximizedGPSFunction( sigma, nu, a, drifters);
optimalAcceleration = fminsearch( errorFunction, log10(0.7e-5), optimset('TolFun',1e-2) )
% a = 3.8915e-06;

% This is a really good looking fit.
% Maximized with a tension on acceleration
nu = 10; sigma_gps = 1.75; sigma = sigma_gps; S = 3;
errorFunction = @(a) StudentTFitMaximizedGPSFunction( sigma, nu, a, drifters);
optimalAcceleration = fminsearch( errorFunction, log10(0.7e-5), optimset('TolFun',1e-2) )
a = 4.7653e-06;

% This was an attempt to repeat the good looking S=3 fit, with S=10.
nu = 10; sigma_gps = 1.75; sigma = sigma_gps; S = 10;
errorFunction = @(a) StudentTFitMaximizedGPSFunction( sigma, nu, a, drifters);
optimalAcceleration = fminsearch( errorFunction, log10(4.7653e-06), optimset('TolFun',1e-2) )
a = 5.16e-06;


% This was run by applying tension on jerk, adn teh results are maybe not
% as good.
nu = 10; sigma_gps = 1.75; sigma = sigma_gps;
errorFunction = @(a) StudentTFitMaximizedGPSFunction( sigma, nu, a, drifters);
optimalAcceleration = fminsearch( errorFunction, log10(1e-9), optimset('TolFun',1e-2) )
a = 5.9113e-09;

% nu = 2.00; sigma = 1.7;
% errorFunction = @(a) StudentTFitScalarFunction( sigma, nu, a, drifters);
% optimalAcceleration = fminsearch( errorFunction, log10(1.0e-5) )
% a = 9.3254e-06;

% nu = 2.00; sigma = 1.25; a = 7.2975e-06;
% errorFunction = @(a) StudentTFitScalarFunction( sigma, nu, a, drifters);
% optimalAcceleration = fminsearch( errorFunction, log10(0.9e-5) )

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
