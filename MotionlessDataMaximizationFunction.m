% params should be [sigma nu]
function lambda = MotionlessDataMaximizationFunction( params, sortedErrors, mu )

sigma = params(1);
nu = params(2);

x = sortedErrors;
n = length(x);
y_data = (1:n)'/n;
y = tcdf((x-mu)/sigma,nu);

D = max(abs(y-y_data));

lambda = (sqrt(n) + 0.12 + 0.11/sqrt(n))*D;
% j=1:25;
% p = sum(2*((-1).^(j-1)).*exp(-2*j.*j*lambda*lambda))