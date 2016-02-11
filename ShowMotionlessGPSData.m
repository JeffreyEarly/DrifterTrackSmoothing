load('sample_data/motionless_gps.mat')

% The GPS was motionless, so its position is the errors
errors = [x;y];


% Now we find the optimal parameters for a student t-distribution
sortedErrors = sort([x;y]);
n = length(sortedErrors);
sigma = std(sortedErrors);
mu = mean(sortedErrors);

errorFunction = @(params) MotionlessDataMaximizationFunction( params, sortedErrors, mu);
[optParams, lambda_t] = fminsearch( errorFunction, [sigma 10], optimset('TolX',1e-1, 'MaxFunEvals', 1e5, 'MaxFunEvals', 1e5) );
sigma_t = optParams(1);
nu = optParams(2);

% Compute the probability this this function is a match
j=1:25;
p_t = sum(2*((-1).^(j-1)).*exp(-2*j.*j*lambda_t*lambda_t))

% Now do the same for the Gaussian distribution
gaussian_cdf = @(z) 0.5*(1 + erf((z-mu)/(sigma*sqrt(2))));  
D = max(abs(gaussian_cdf(sortedErrors) - ((1:n)'/n)));
lambda_g = (sqrt(n) + 0.12 + 0.11/sqrt(n))*D;
p_g = sum(2*((-1).^(j-1)).*exp(-2*j.*j*lambda_g*lambda_g))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Make some plots!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(2,1,1)
histogram(errors,100)

subplot(2,1,2)
plot(sortedErrors,(1:n)'/n, 'k', 'LineWidth', 2)
hold on
z = linspace(-30,30,100)';
plot(z, gaussian_cdf(z), 'b', 'LineWidth', 1)
plot(z, tcdf((z-mu)/sigma_t,nu), 'g', 'LineWidth', 1)
legend('data','Gaussian','Student t')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Okay, now analyze these separately
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sortedXErrors = sort(x);
sortedYErrors = sort(y);
n = length(sortedXErrors);
mu_x = mean(sortedErrors);
mu_y = mean(sortedErrors);

errorFunction = @(params) MotionlessDataMaximizationFunction( params, sortedXErrors, mu_x);
[optXParams, lambdaX_t] = fminsearch( errorFunction, [std(x) nu], optimset('TolX',1e-1, 'MaxFunEvals', 1e5, 'MaxFunEvals', 1e5) )

errorFunction = @(params) MotionlessDataMaximizationFunction( params, sortedYErrors, mu_y);
[optYParams, lambdaY_t] = fminsearch( errorFunction, [std(y) nu], optimset('TolX',1e-1, 'MaxFunEvals', 1e5, 'MaxFunEvals', 1e5) )


figure
subplot(2,2,1)
histogram(x,100)

subplot(2,2,2)
histogram(y,100)

subplot(2,2,3)
plot(sortedXErrors,(1:n)'/n, 'k', 'LineWidth', 2)
hold on
z = linspace(-30,30,100)';
plot(z, tcdf((z-mu_x)/optXParams(1),optXParams(2)), 'g', 'LineWidth', 2)
legend('data','Student t')

subplot(2,2,4)
plot(sortedYErrors,(1:n)'/n, 'k', 'LineWidth', 2)
hold on
z = linspace(-30,30,100)';
plot(z, tcdf((z-mu_y)/optYParams(1),optYParams(2)), 'g', 'LineWidth', 2)
legend('data','Student t')