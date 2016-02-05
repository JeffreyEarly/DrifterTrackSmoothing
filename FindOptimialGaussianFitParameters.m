drifters = load('sample_data/projected_ungridded_rho1_drifters.mat');

iDrifter = 6;

x = drifters.x{iDrifter};
y = drifters.y{iDrifter};
t = drifters.t{iDrifter};

S = 3; % order of the spline
K = S+1;

errorFunction = @(x0) GaussianFitScalarFunctionBivariate(x0,t,x,y,S,'total');
optimalParameters = fminsearch( errorFunction, [150, 1.0e-5] )
% sigma = 80;
% a = 1.1826e-04;

errorFunction = @(x0) GaussianFitScalarFunctionBivariate(x0,t,x,y,S,'acceleration');
optimalParameters = fminsearch( errorFunction, [150, log10(1.0e-4)] ) % [77.8331   -3.8579]
optimalParameters = fminsearch( errorFunction, [150, log10(1.0e-5)] ) % [131.8940   -5.6154]
optimalParameters = fminsearch( errorFunction, [150, log10(1.0e-6)] ) % [159.8404   -5.7565];
% a = 1.075e-5;
% sigma = 66;

% We choose a value of acceleration that is consistent with u_rms*f0
u_rms = 0.2;
a = u_rms*drifters.f0;

% Find the optimal fit to acceleration, varying sigma
a = 1.5e-5;
errorFunction = @(sigma) GaussianFitScalarFunction(sigma, a,t,x,y,S,'acceleration');
% optimalSigma = fminsearch( errorFunction, 150 )
% optimalSigma = 65

% Find the optimal fit to position, varying sigma
a = 1.5e-5;
errorFunction = @(sigma) GaussianFitScalarFunction(sigma, a,t,x,y,S,'position');
% optimalSigma = fminsearch( errorFunction, 150 )
% optimalSigma = 195

% Find the optimal fit to position + acceleration, varying sigma
a = 1.5e-5;
errorFunction = @(sigma) GaussianFitScalarFunction(sigma, a,t,x,y,S,'total');
% optimalSigma = fminsearch( errorFunction, 150 )
% optimalSigma = 195
% We find that the total error follows the position error. Thus, the error
% is dominanted by positions errors.


% % Find the optimal fit to sigma, varying sigma
sigma = 200;
% errorFunction = @(a) GaussianFitScalarFunction(sigma, a,t,x,y,S,'acceleration');
% optimalAccleration = fminsearch( errorFunction, 1e-5 )

errorFunction = @(a) GaussianFitScalarFunction(sigma, a,t,x,y,S,'position');
optimalAccleration = fminsearch( errorFunction, 1e-5 )

errorFunction = @(a) GaussianFitScalarFunction(sigma, a,t,x,y,S,'total');
optimalAccleration = fminsearch( errorFunction, 1e-5 )