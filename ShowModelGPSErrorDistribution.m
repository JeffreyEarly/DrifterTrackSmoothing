sigma_gps = 1.5; nu = 6; % 95 percent is at 4.75 m, 99.99 percent is at 15.5 meters
sigma_gps = 1.75; nu = 5; % 95 percent is at 5.88 m, 99.99 percent is at 22.66 meters
sigma_gps = 1.0; nu = 6; % 95 percent is at 3.15 m, 99.99 percent is at 10.3 meters
sigma_gps = 1.25; nu = 10; % 95 percent is at 3.54 m, 99.99 percent is at 5.58 meters

sigma_gps = 8.0; nu = 5.5; % 95 percent is at 3.54 m, 99.99 percent is at 5.58 meters

pdf = @(x) (x/sigma_gps^2).*exp(-x.*x/(2*sigma_gps^2));
cdf = @(x) 1-exp(-x.*x/(1.75^2));

x = linspace(0,20,1000)';
cdf1 = cdf(x); 
v95 = x(find( cdf1 < 0.95, 1, 'last'))
v9999 = x(find( cdf1 < 0.9995, 1, 'last'))

x = linspace(0,20,1000)';
figure
plot(x, pdf(x));
vlines([v95, v9999],'g--')

sigma = sigma_gps;

pdf = @(z) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2));
cdf = @(z) 1/2 + x.* gamma((nu+1)/2) * HyperGeometric2F1([1/2; (nu+1)/2],3/2,-z/(sigma*sigma*nu))/( sqrt(pi*nu)*gamma(nu/2));

m = 100;
x = linspace(-m,m,2000)';
y = linspace(-m,m,2000)';
[X,Y] = meshgrid(x,y);
rho = sqrt( X.*X + Y.*Y);
dR = x(2)-x(1);

pdf2d = pdf(X).*pdf(Y);
pdf2d_norm = pdf2d * dR;

rMag = 0:dR:max(max(rho));
indices = cell( length(rMag), 1);
for i = 1:length(rMag)
	indices{i} = find( rho >= rMag(i)-dR/2 & rho < rMag(i)+dR/2 );
end

% Now sum up the energy
pdf1d = zeros( size(rMag) );
for i = 1:length(rMag)
	pdf1d(i) = sum( pdf2d_norm( indices{i} ) );
end

cdf1d = cumtrapz(rMag,pdf1d);
v95 = rMag(find( cdf1d < 0.95, 1, 'last'))
v9999 = rMag(find( cdf1d < 0.9995, 1, 'last'))

figure
plot(rMag, pdf1d)
vlines([v95, v9999],'g--')