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

iDrifter = 8;
SplineFactor = 3.0; % Number of data points for each spline
sigma_gps = 10; % error in meters
S = 5; % order of the spline
u_rms = 0.15; % assumed rms velocity of the solution
T_decorrelation = 0; %2.5*60*60; % forcing decorrelation time

% Note that big decorrelation times are causing the the solvability issue.
% for some reason the the x3 points are not collocating.

lat0 = drifters.lat0;
x = drifters.x{iDrifter};
y = drifters.y{iDrifter};
t = drifters.t{iDrifter};

error_y = ones(size(y))*sigma_gps; % This is the error, in kilometers, of the GPS fits.
N = length(y);
M = floor(N/SplineFactor); % Number of splines

dx = ones(size(x))*sigma_gps;
dy = ones(size(y))*sigma_gps;


% [mx,my,Cmx,Cmy,A,V] = drifter_fit(t,x,y,dx,dy,W,M,1000,0, @(z)(z./(1+0.5*z.*z)));
[mx,my,Cmx,Cmy,A,V] = drifter_fit_forced(t,x,y,dx,dy,T_decorrelation,M,S,1000,0,@(z)(z));
% [mx,my,Cmx,Cmy,A,V] = drifter_fit_generalized(t,x,y,dx,dy,M,T_decorrelation,1000,lat0, @(z)(z./(1+0.5*z.*z)));
x1 = A*mx;
y1 = A*my;
u1 = V*mx;
v1 = V*my;


[mx,my,Cmx,Cmy,X2,V] = drifter_fit_forced(t,x,y,dx,dy,T_decorrelation,M,S,u_rms,0, @(z)(z./(1+0.5*z.*z)));
% [mx,my,Cmx,Cmy,A,V] = drifter_fit_generalized(t,x,y,dx,dy,M,T_decorrelation,u_rms,0, @(z)(z./(1+0.5*z.*z)));
x2 = X2*mx;
y2 = X2*my;
u2 = V*mx;
v2 = V*my;



[mx,my,Cmx,Cmy,X2,V] = drifter_fit_forced(t,x,y,dx,dy,T_decorrelation,M,S,u_rms,lat0, @(z)(z./(1+0.5*z.*z)));
% [mx,my,Cmx,Cmy,A,V] = drifter_fit_generalized(t,x,y,dx,dy,M,T_decorrelation,u_rms,lat0, @(z)(z./(1+0.5*z.*z)));
x3 = X2*mx;
y3 = X2*my;
u3 = V*mx;
v3 = V*my;

% Now we create the basis at the desired collocation points.
t = (0:30*60:drifters.maxExperimentLength)';
Nt = length(t);
if S == 3
    addpath('./cubic_splines');
    spline = @(t) cspline(t);
    spline_t = @(t) cspline_t(t);
    spline_tt = @(t) cspline_tt(t);
    spline_ttt = @(t) cspline_ttt(t);
elseif S == 5
    addpath('./quintic_splines');
    spline = @(t) qspline(t);
    spline_t = @(t) qspline_t(t);
    spline_tt = @(t) qspline_tt(t);
    spline_ttt = @(t) qspline_ttt(t);
else
    disp('Whoops! I only know how to deal with splines of order 3 and 5.')
    return;
end
M_norm = M+2*floor(S/2);
t_knot = (drifters.t{iDrifter}(end)-drifters.t{iDrifter}(1))/(M_norm-S);
X = zeros(Nt,M_norm);
V = zeros(Nt,M_norm);
A = zeros(Nt,M_norm);
J = zeros(Nt,M_norm);
for i=1:Nt
    for j=1:M_norm
        t_norm=(t(i)-drifters.t{iDrifter}(1))/t_knot - (j - 1 - floor(S/2));
        X(i,j)=spline(t_norm);
        V(i,j)=spline_t(t_norm);
        A(i,j)=spline_tt(t_norm);
        J(i,j)=spline_ttt(t_norm);
    end
end
V = V/t_knot;
A = A/t_knot^2;
J = J/t_knot^3;

% Compute the Coriolis parameter
Omega = 2*pi/86164;
f0 = 2*Omega*sin(lat0*pi/180);

x4 = X*mx;
y4 = X*my;
u4 = V*mx;
v4 = V*my;
ax4 = A*mx;
ay4 = A*my;

f_x = ax4 + f0*v4;
f_y = ay4 - f0*u4;
ft_x = J*mx + f0*A*my;
ft_y = J*my - f0*A*mx;

figure
subplot(2,1,1)
plot(t/3600,[f_x,f_y])
subplot(2,1,2)
plot(t/3600,[ft_x,ft_y])

figure
subplot(2,1,1)
plot(t/3600,[A*mx,A*my])
subplot(2,1,2)
plot(t/3600,[J*mx,J*my])

figure
subplot(2,2,[1 3])
s = 1/1000;
plot(s*x,s*y), hold on
plot(s*x1,s*y1,'g')
plot(s*x2,s*y2,'r')
plot(s*x3,s*y3,'k')
plot(s*x4,s*y4,'c')
scatter(s*x,s*y,5)
% xlim([-500 7500])
% ylim([2000 10000])
xlabel('x (km)')
ylabel('y (km)')
title('Fit with Cauchy errors, no tension')

subplot(2,2,2)
plot(drifters.t{iDrifter}/3600,s*x), hold on
plot(drifters.t{iDrifter}/3600,s*x1,'g')
plot(drifters.t{iDrifter}/3600,s*x2,'r')
plot(drifters.t{iDrifter}/3600,s*x3,'k')
plot(t/3600,s*x4,'c')
scatter(drifters.t{iDrifter}/3600,s*x,5)
% xlim([8.295e6 1.01e7])
xlabel('t (hours)')
ylabel('x (km)')

subplot(2,2,4)
plot(drifters.t{iDrifter}/3600,s*y), hold on
plot(drifters.t{iDrifter}/3600,s*y1,'g')
plot(drifters.t{iDrifter}/3600,s*y2,'r')
plot(drifters.t{iDrifter}/3600,s*y3,'k')
plot(t/3600,s*y4,'c')
scatter(drifters.t{iDrifter}/3600,s*y,5)
% xlim([8.295e6 1.01e7])
xlabel('t (hours)')
ylabel('y (km)')

% output = sprintf('%s/center_of_mass.eps', drifters.figuresFolder);
% print('-depsc2', output)

xi = (-0500:0.01:0500)';

dx = x1 - x;
dy = y1 - y;
ds = sqrt(dx.*dx+dy.*dy);
sigma_x = std(dx,1);
sigma_y = std(dy,1);
sigma_s = std(ds,1);

nsigma = 2;

figure
subplot(1,3,1)
denhist(X2*mx - x, 500,'b'); hold on
plot(xi,exp(-(xi.*xi)/(2*sigma_gps*sigma_gps))/(sigma_gps*sqrt(2*pi)),'LineWidth',2)
plot(xi,1./(pi*sigma_gps*(1+(xi.*xi)/(sigma_gps*sigma_gps))),'LineWidth',2)
vlines([-nsigma*sigma_x nsigma*sigma_x],'r--')
xlim([-250 250])
xlabel('model x-error (meters)')

subplot(1,3,2)
denhist(X2*my - y, 500,'b'); hold on
plot(xi,exp(-(xi.*xi)/(2*sigma_gps*sigma_gps))/(sigma_gps*sqrt(2*pi)),'LineWidth',2)
plot(xi,1./(pi*sigma_gps*(1+(xi.*xi)/(sigma_gps*sigma_gps))),'LineWidth',2)
vlines([-nsigma*sigma_y nsigma*sigma_y],'r--')
xlim([-250 250])
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

