%clear
drifters = load('sample_data/projected_ungridded_rho1_drifters.mat');

% These drifters are already projected, so we don't need this.
%
% t = date{iDrifter}-date{iDrifter}(1);
% lat0 = min(lat{iDrifter}) + (max(lat{iDrifter})-min(lat{iDrifter}))/2;
% lon0 = min(lon{iDrifter}) + (max(lon{iDrifter})-min(lon{iDrifter}))/2;
% [x,y] = latlon2xy( lat{iDrifter}, lon{iDrifter}, lat0, lon0 );
% t = t*24*60*60; % convert from days to seconds
% x = x*1000; % convert from km to meters.
% y = y*1000;

iDrifter = 1;
sigma_gps = 9; % error in meters
nu = 2.015;
S = 3; % order of the spline
u_rms = 0.15; % assumed rms velocity of the solution

% Note that big decorrelation times are causing the the solvability issue.
% for some reason the the x3 points are not collocating.

lat0 = drifters.lat0;
x = drifters.x{iDrifter};
y = drifters.y{iDrifter};
t = drifters.t{iDrifter};

dx = ones(size(x))*sigma_gps;
dy = ones(size(y))*sigma_gps;

% t_knot = t(1:3:end);
% t_knot = [t(1); t(4:3:end-4); t(end)];
t_knot = t;

% Gaussian and t-distribution
p_g = @(z,sigma) exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi));
p_t = @(z,sigma) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2));

% Gaussian and t-distribution weighting functions.
w_g = @(z)(sigma_gps*sigma_gps);
w_t = @(z)((nu/(nu+1))*sigma_gps^2*(1+z.^2/(nu*sigma_gps^2)));

% [mx,my,Cmx,Cmy,A,V] = drifter_fit(t,x,y,dx,dy,W,M,1000,0, @(z)(z./(1+0.5*z.*z)));
[m_x,m_y,Cm_x,Cm_y,B,Bq] = drifter_fit_bspline(t,x,y,dx,dy,S,t_knot,[0,1.5e12],w_g);
X = squeeze(B(:,:,1));
Xq = squeeze(Bq(:,:,1));
Vq = squeeze(Bq(:,:,2));

x1 = Xq*m_x;
y1 = Xq*m_y;
u1 = Vq*m_x;
v1 = Vq*m_y;
tq = linspace(drifters.t{iDrifter}(1), drifters.t{iDrifter}(end),size(Xq,1));

[m_x,m_y,Cm_x,Cm_y,B,Bq] = drifter_fit_bspline(t,x,y,dx,dy,S,t_knot,[0,6.25e10],w_t);
X = squeeze(B(:,:,1));
Xq = squeeze(Bq(:,:,1));
Vq = squeeze(Bq(:,:,2));
Aq = squeeze(Bq(:,:,3));
Jq = squeeze(Bq(:,:,4));
x2 = Xq*m_x;
y2 = Xq*m_y;
u2 = Vq*m_x;
v2 = Vq*m_y;
jx2 = Jq*m_x;
jy2 = Jq*m_y;

Smoothness = 2.5;

M = floor(length(t)/Smoothness);
j2 = sqrt(jx2.*jx2+jy2.*jy2);
xi = cumsum(abs(j2).^(1/(S+1)))*(tq(2)-tq(1));
xi_interp = linspace(xi(1),xi(end),M);
t_knot = interp1(xi,tq,xi_interp,'linear')';

[m_x,m_y,Cm_x,Cm_y,B,Bq] = drifter_fit_bspline(t,x,y,dx,dy,S,t_knot,10,w_t);
x3 = Xq*m_x;
y3 = Xq*m_y;
u3 = Vq*m_x;
v3 = Vq*m_y;
jx3 = Jq*m_x;
jy3 = Jq*m_y;

M = floor(length(t)/Smoothness);
j3 = sqrt(jx3.*jx3+jy3.*jy3);
xi = cumsum(abs(j3).^(1/(S+1)))*(tq(2)-tq(1));
xi_interp = linspace(xi(1),xi(end),M);
t_knot3 = interp1(xi,tq,xi_interp,'linear')';

[m_x,m_y,Cm_x,Cm_y,B,Bq] = drifter_fit_bspline(t,x,y,dx,dy,S,t_knot3,10,w_t);
x3 = Xq*m_x;
y3 = Xq*m_y;
u3 = Vq*m_x;
v3 = Vq*m_y;
jx3 = Jq*m_x;
jy3 = Jq*m_y;

M = floor(length(t)/Smoothness);
j3 = sqrt(jx3.*jx3+jy3.*jy3);
xi = cumsum(abs(j3).^(1/(S+1)))*(tq(2)-tq(1));
xi_interp = linspace(xi(1),xi(end),M);
t_knot3 = interp1(xi,tq,xi_interp,'linear')';

[m_x,m_y,Cm_x,Cm_y,X,V,A,J,Xq,Vq,Aq,Jq] = drifter_fit_bspline(t,x,y,dx,dy,S,t_knot3,10,w_t);
x3 = Xq*m_x;
y3 = Xq*m_y;
u3 = Vq*m_x;
v3 = Vq*m_y;
jx3 = Jq*m_x;
jy3 = Jq*m_y;

figure
subplot(2,2,[1 3])
s = 1/1000;
plot(s*x,s*y), hold on
plot(s*x1,s*y1,'g')
plot(s*x2,s*y2,'r')
plot(s*x3,s*y3,'k')
% plot(s*x4,s*y4,'c')
scatter(s*x,s*y,5)
% scatter(s*interp1(tq,x3,t_knot),s*interp1(tq,y3,t_knot),5)
% xlim([-500 7500])
% ylim([2000 10000])
xlabel('x (km)')
ylabel('y (km)')
title('Fit with Cauchy errors, no tension')

subplot(2,2,2)
plot(drifters.t{iDrifter}/3600,s*x), hold on
plot(tq/3600,s*x1,'g')
plot(tq/3600,s*x2,'r')
plot(tq/3600,s*x3,'k')
% plot(t/3600,s*x4,'c')
scatter(drifters.t{iDrifter}/3600,s*x,5)
% scatter(t_knot/3600,s*interp1(tq,x3,t_knot),5)
% xlim([8.295e6 1.01e7])
xlabel('t (hours)')
ylabel('x (km)')

subplot(2,2,4)
plot(drifters.t{iDrifter}/3600,s*y), hold on
plot(tq/3600,s*y1,'g')
plot(tq/3600,s*y2,'r')
plot(tq/3600,s*y3,'k')
% plot(t/3600,s*y4,'c')
scatter(drifters.t{iDrifter}/3600,s*y,5)
% scatter(t_knot/3600,s*interp1(tq,y3,t_knot),5)
% xlim([8.295e6 1.01e7])
xlabel('t (hours)')
ylabel('y (km)')

xi = (-0500:0.01:0500)';

dx = X*m_x - x;
dy = X*m_y - y;
ds = sqrt(dx.*dx+dy.*dy);
sigma_x = std(dx,1);
sigma_y = std(dy,1);
sigma_s = std(ds,1);

ACx = Autocorrelation(dx,30);
ACy = Autocorrelation(dy,30);

figure, plot(ACx(2:end))



nsigma = 2;

figure
subplot(1,3,1)
denhist(X*m_x - x, 500,'b'); hold on
plot(xi,exp(-(xi.*xi)/(2*sigma_gps*sigma_gps))/(sigma_gps*sqrt(2*pi)),'LineWidth',2)
plot(xi,1./(pi*sigma_gps*(1+(xi.*xi)/(sigma_gps*sigma_gps))),'LineWidth',2)
vlines([-nsigma*sigma_x nsigma*sigma_x],'r--')
% xlim([-250 250])
xlim([-50 50])
xlabel('model x-error (meters)')

subplot(1,3,2)
denhist(X*m_y - y, 500,'b'); hold on
plot(xi,exp(-(xi.*xi)/(2*sigma_gps*sigma_gps))/(sigma_gps*sqrt(2*pi)),'LineWidth',2)
plot(xi,1./(pi*sigma_gps*(1+(xi.*xi)/(sigma_gps*sigma_gps))),'LineWidth',2)
vlines([-nsigma*sigma_y nsigma*sigma_y],'r--')
% xlim([-250 250])
xlim([-50 50])
xlabel('model y-error (meters)')

subplot(1,3,3)
denhist(ds,500);
xlabel('model distance error (meters)')

figure
subplot(4,2,1)
hist(diff(x)./diff(drifters.t{iDrifter}),50)
title(sprintf('u raw, rms=%f',sqrt(mean((diff(x)./diff(drifters.t{iDrifter})).^2))))
subplot(4,2,2)
hist(diff(y)./diff(drifters.t{iDrifter}),50)
title(sprintf('v raw, rms=%f',sqrt(mean((diff(y)./diff(drifters.t{iDrifter})).^2))))
subplot(4,2,3)
hist(u1,50)
title(sprintf('u fit, rms=%f',sqrt(mean(u1.^2))))
subplot(4,2,4)
hist(v1,50)
title(sprintf('v fit, rms=%f',sqrt(mean(v1.^2))))
subplot(4,2,5)
hist(u2,50)
title(sprintf('u fit tension, rms=%f',sqrt(mean(u2.^2))))
subplot(4,2,6)
hist(v2,50)
title(sprintf('v fit tension, rms=%f',sqrt(mean(v2.^2))))
subplot(4,2,7)
hist(u3,50)
title(sprintf('u fit tension+inertial, rms=%f',sqrt(mean(u3.^2))))
subplot(4,2,8)
hist(v3,50)
title(sprintf('v fit tension+inertial, rms=%f',sqrt(mean(v3.^2))))

