%clear
load('All2011Drifters.mat')
addpath('../support/splines/')

preferredIndices = rho1DrifterIndices;

% The last drifter has a shorter time series.
% preferredIndices(end)=[];
% preferredIndices(1)=[];

lat = {lat{preferredIndices}};
lon = {lon{preferredIndices}};
date = {date{preferredIndices}};

numDrifters = length(date);

%
%	Determine the first deployment date
%
lastDeployment = date{1}(1);
firstRetrieval = date{1}(end);
for i=1:numDrifters
	if ( date{i}(1) > lastDeployment )
		lastDeployment = date{i}(1);
	end
	if ( date{i}(end) < firstRetrieval )
		firstRetrieval = date{i}(end);
	end
end

%
% Filter the raw data by removing points with spikes in curvature
%

iDrifter = 2;
SplineFactor = 1.0;
gamma = 0.0000; % Tension
sigma_gps = 15; % error in meters

t = date{iDrifter}-date{iDrifter}(1);
lat0 = min(lat{iDrifter}) + (max(lat{iDrifter})-min(lat{iDrifter}))/2;
lon0 = min(lon{iDrifter}) + (max(lon{iDrifter})-min(lon{iDrifter}))/2;
[x,y] = latlon2xy( lat{iDrifter}, lon{iDrifter}, lat0, lon0 );

t = t*24*60*60; % convert from days to seconds
x = x*1000; % convert from km to meters.
y = y*1000;

error_y = ones(size(y))*sigma_gps; % This is the error, in kilometers, of the GPS fits.
N = length(y);
M = floor(N/SplineFactor); % Number of splines

dx = ones(size(x))*sigma_gps;
dy = ones(size(y))*sigma_gps;

u_rms = 0.09;
T_decorrelation = 5*60*60;
% T_decorrelation = 1;

W=zeros(N,N);
for i=1:N
    for j=1:N
        W(i,j)=abs(t(i)-t(j));
    end
end
% T_g2 = -T_decorrelation*T_decorrelation/log(0.01);
% W=exp(-W.*W/T_g2);

T_e = -T_decorrelation/log(0.01);
W=exp(-W/T_e);

% W2 = W;
% W2(find(W2<0.01))=0.0;
% return;
% W=eye(N);

[mx,my,Cmx,Cmy,A,V] = drifter_fit(t,x,y,dx,dy,W,M,u_rms,lat0, @(z)(z./(1+0.5*z.*z)));
% [mx,my,Cmx,Cmy,A,V] = drifter_fit_generalized(t,x,y,dx,dy,M,T_decorrelation,u_rms,lat0, @(z)(z./(1+0.5*z.*z)));
x3 = A*mx;
y3 = A*my;
u3 = V*mx;
v3 = V*my;

% [mx,my,Cmx,Cmy,A,V] = drifter_fit(t,x,y,dx,dy,W,M,1000,0, @(z)(z./(1+0.5*z.*z)));
[mx,my,Cmx,Cmy,A,V] = drifter_fit(t,x,y,dx,dy,W,M,1000,0,@(z)(z));
% [mx,my,Cmx,Cmy,A,V] = drifter_fit_generalized(t,x,y,dx,dy,M,T_decorrelation,1000,lat0, @(z)(z./(1+0.5*z.*z)));
x1 = A*mx;
y1 = A*my;
u1 = V*mx;
v1 = V*my;

[mx,my,Cmx,Cmy,A,V] = drifter_fit(t,x,y,dx,dy,W,M,u_rms,0, @(z)(z./(1+0.5*z.*z)));
% [mx,my,Cmx,Cmy,A,V] = drifter_fit_generalized(t,x,y,dx,dy,M,T_decorrelation,u_rms,0, @(z)(z./(1+0.5*z.*z)));
x2 = A*mx;
y2 = A*my;
u2 = V*mx;
v2 = V*my;



figure
subplot(2,2,[1 3])
plot(x,y)
hold on
plot(x1,y1,'g')
plot(x2,y2,'r')
plot(x3,y3,'k')
scatter(x,y,5)
% xlim([-500 7500])
% ylim([2000 10000])
xlabel('x (meters)')
ylabel('y (meters)')
title('Fit with Cauchy errors, no tension')

subplot(2,2,2)
plot(t*24,x), hold on
plot(t*24,x1,'g')
plot(t*24,x2,'r')
plot(t*24,x3,'k')
scatter(t*24,x,5)
% xlim([8.295e6 1.01e7])
xlabel('t (seconds)')
ylabel('x (meters)')

subplot(2,2,4)
plot(t*24,y), hold on
plot(t*24,y1,'g')
plot(t*24,y2,'r')
plot(t*24,y3,'k')
scatter(t*24,y,5)
% xlim([8.295e6 1.01e7])
xlabel('t (seconds)')
ylabel('y (meters)')

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
denhist(A*mx - x, 500,'b'); hold on
plot(xi,exp(-(xi.*xi)/(2*sigma_gps*sigma_gps))/(sigma_gps*sqrt(2*pi)),'LineWidth',2)
plot(xi,1./(pi*sigma_gps*(1+(xi.*xi)/(sigma_gps*sigma_gps))),'LineWidth',2)
vlines([-nsigma*sigma_x nsigma*sigma_x],'r--')
xlim([-250 250])
xlabel('model x-error (meters)')

subplot(1,3,2)
denhist(A*my - y, 500,'b'); hold on
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
hist(diff(x)./diff(t),50)
title(sprintf('u raw, rms=%f',sqrt(mean((diff(x)./diff(t)).^2))))
subplot(4,2,2)
hist(diff(y)./diff(t),50)
title(sprintf('v raw, rms=%f',sqrt(mean((diff(y)./diff(t)).^2))))
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

