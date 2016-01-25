% close all, clear

% c = []; d=[]; f=[]; g=[];



% Create a simple signal
dt = 1;
t=(0:dt:100)'; % seconds
a_true = 0;
u_true = 1; % meters/second

path = @(t) a_true*(t.*t) + u_true.*t;
speed = @(t) 2*a_true.*t + u_true*ones(size(t));

x_true = path(t);

sigma = 4; % meters

% Create some Gaussian noise
epsilon = sigma*randn(size(x_true));

% Contaminate the signal
x = x_true + epsilon;


p = @(z) exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi));
w = @(z)(sigma*sigma);

S=3;
t_knot = t;
N = length(t);

v = 1e-3;
gamma = (N)./(v*v*(t(end)-t(1)));
[m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = drifter_fit_bspline(t,x,x,ones(size(x))*sigma,ones(size(x))*sigma,S,t_knot,gamma,w);
x_fit_low_tension = squeeze(Bq(:,:,1))*m_x;

% range = 10.^linspace(-5,-1,10);

% for i=1:length(range)
v = 1*sigma/dt;
a = 0.1*sqrt(2)*sigma/(dt*dt);
gamma_v = (N)./(v*v*(t(end)-t(1)));
gamma_a = (N)./(a*a*(t(end)-t(1)));
gamma = [gamma_v;0*gamma_a];
[m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = drifter_fit_bspline(t,x,x,ones(size(x))*sigma,ones(size(x))*sigma,S,t_knot,gamma,w);
x_fit = squeeze(Bq(:,:,1))*m_x;
v_fit = squeeze(Bq(:,:,2))*m_x;
a_fit = squeeze(Bq(:,:,3))*m_x;
j_fit = squeeze(Bq(:,:,4))*m_x;
x_error = x - squeeze(B(:,:,1))*m_x;

% c(end+1)=a;
% d(end+1)=std(v_fit);
% f(end+1)=std(a_fit);
% g(end+1)=std(j_fit);

% end

mean_x_error = sqrt(mean((path(tq) - x_fit).^2));
mean_v_error = sqrt(mean((speed(tq) - v_fit).^2));

figure
plot(t,x_true), hold on
plot(tq,x_fit_low_tension,'g')
plot(tq,x_fit,'b')
scatter(t,x,5)
title(sprintf('position rms error=%.2f meters, velocity rms error=%.2f m/s',mean_x_error, mean_v_error))

figure
hist(x_error)
title(sprintf('std=%f',std(x_error)))
