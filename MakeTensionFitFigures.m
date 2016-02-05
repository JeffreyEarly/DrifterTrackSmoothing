% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

drifters = load('sample_data/projected_ungridded_rho1_drifters.mat');

 distribution = 'student-t';
 distribution = 'gaussian';
 
 iDrifter = 6;


S = 2; % order of the spline
K = S+1;

if strcmp(distribution,'gaussian')
    sigma = 150; % error in meters
    
    p = @(z) exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi));
    w = @(z)(sigma*sigma);
    
    a = 1.4e-5;
    j = 7e-9;
    gamma = [0; 1/a^2; 1/j^2];
    
elseif strcmp(distribution,'student-t')
    % nu = 2.005; sigma = 8;
    nu = 2.019; sigma = 15;
    %var = 158.^2; sigma = 20; nu = 2*var/(var-sigma*sigma);
    
    a = 1.4e-6;
    j = 7e-9;
    gamma = [0; 1/a^2; 1/j^2];
     
    p = @(z) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2));
    w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));
    rho = @(z) -log(sqrt(2)*gamma((nu+1)/2)./(sqrt(nu)*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2)));
else
    disp('Not a valid distribution')
    return;
end

x = drifters.x{iDrifter};
y = drifters.y{iDrifter};
t = drifters.t{iDrifter};

dx = ones(size(x))*sigma;
dy = ones(size(y))*sigma;

N = length(t);

[m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = bspline_bivariate_fit_with_tension(t,x,y,dx,dy,S,gamma, w);


Xq = squeeze(Bq(:,:,1));
x_fit = Xq*m_x;
y_fit = Xq*m_y;

figure
subplot(2,2,[1 3])
s = 1/1000;
plot(s*x,s*y), hold on
plot(s*x_fit,s*y_fit,'g')
scatter(s*x,s*y,5)
xlabel('x (km)')
ylabel('y (km)')

subplot(2,2,2)
plot(drifters.t{iDrifter}/3600,s*x), hold on
plot(tq/3600,s*x_fit,'g')
scatter(drifters.t{iDrifter}/3600,s*x,5)
xlabel('t (hours)')
ylabel('x (km)')

subplot(2,2,4)
plot(drifters.t{iDrifter}/3600,s*y), hold on
plot(tq/3600,s*y_fit,'g')
scatter(drifters.t{iDrifter}/3600,s*y,5)
xlabel('t (hours)')
ylabel('y (km)')