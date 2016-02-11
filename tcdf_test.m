
sigma = 1.75;
nu = 5;

x = linspace(-200,200,1000)';

pdf = @(z) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2));
cdf = @(z) tcdf(z,nu);

figure
plot(x,cumtrapz(x,pdf(x)))
hold on
plot(x,cdf(x/sigma))
xlim([-20 20])