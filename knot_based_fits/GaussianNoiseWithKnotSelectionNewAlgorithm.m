% close all, clear

% c = []; d=[]; f=[]; g=[];

addpath('../support')


rng(3)

%%%%%%%%%%%%%%%%%%%%%%%%%%
% quadratic
%%%%%%%%%%%%%%%%%%%%%%%%%%
dt_v = 1;
t=(0:dt_v:100)'; % seconds
a_true = .1;
u_true = 0; % meters/second

% path = @(t) a_true*(t.*t) + u_true.*t;
% speed = @(t) 2*a_true.*t + u_true*ones(size(t));
% acceleration = @(t) 2*a_true*ones(size(t));

%%%%%%%%%%%%%%%%%%%%%%%%%%
% sech^2
%%%%%%%%%%%%%%%%%%%%%%%%%%
u_true = 500;
u_width = 10;
t_0 = 50;
path = @(t) u_true*sech((t-t_0)/u_width).^2;
speed = @(t) -2*(u_true/u_width).*tanh((t-t_0)/u_width).*sech((t-t_0)/u_width).^2;
acceleration = @(t) (u_true/u_width^2)*(4 * (tanh((t-t_0)/u_width)).^2 .* (sech((t-t_0)/u_width)).^2 - 2*(sech((t-t_0)/u_width)).^4);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Just noise
%%%%%%%%%%%%%%%%%%%%%%%%%%
% path = @(t) zeros(size(t));
% speed = @(t) zeros(size(t));
% acceleration = @(t) zeros(size(t));

x_true = path(t);
v_true = speed(t);
a_true = acceleration(t);

sigma = 50; % meters

% Create some Gaussian noise
epsilon = sigma*randn(size(x_true));

% Contaminate the signal
x = x_true + epsilon;

N = length(t);

p = @(z) exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi));
w = @(z)(sigma*sigma);


z_threshold = 2.4;

Sigma = sigma*ones(size(t));
S = 0;
[t_knot0, S, constraints] = FindStatisticallySignificantKnotRegions(t,x,Sigma,z_threshold,w, S);
[m_x0,Cm_x0,B0] = bspline_fit_no_tension_constrain(t,x,Sigma,S,t_knot0,w, constraints);

tq0 = linspace(t(1),t(end),10*N)';
Bq0 = bspline(tq0,t_knot0,S+1);
tm0 = (t_knot0(1:end-1) + t_knot0(2:end))/2;
Bm0 = bspline(tm0,t_knot0,S+1);
val = squeeze(Bm0(:,:,1))*m_x0;
Ex = squeeze(Bm0(:,:,1))*Cm_x0*squeeze(Bm0(:,:,1)).';
% A == B
% A = diag(Ex);
% B = group0.error;

x_fit0 = squeeze(Bq0(:,:,1))*m_x0;
x_fit0_error = sqrt(diag(squeeze(Bq0(:,:,1))*Cm_x0*squeeze(Bq0(:,:,1)).'));
x_error0 = x - squeeze(B0(:,:,1))*m_x0;
mean_x_error0 = sqrt(mean((path(tq0) - x_fit0).^2));
mean_v_error0 = sqrt(mean((speed(tq0)).^2));
mean_a_error0 = sqrt(mean((acceleration(tq0)).^2));

Q_error0 = (mean_x_error0 / sqrt(mean((path(tq0)).^2)) - 1);

% fprintf('This is how much the error decreases as data points are added.\n')
% S = 2;
% K = S+1;
% % What are the returned errors
% t_knot1 = NaturalKnotsForSpline(t,K);
% loops = 5;
% v = zeros(10,loops);
% for i=1:loops
% [m_x1,Cm_x1,B1] = bspline_fit_no_tension_constrain(t,x,ones(size(x))*sigma,S,t_knot1,w);
% tm1 = (t_knot1(1:end-1) + t_knot1(2:end))/2;
% Bm1 = bspline(tm1,t_knot1,S+1);
% val = squeeze(Bm1(:,:,K))*m_x1;
% Ex = squeeze(Bm1(:,:,K))*Cm_x1*squeeze(Bm1(:,:,K)).';
% tmp = diag(Ex);
% % tmp(1:10)/(sigma*sigma)
% v(:,i) = sigma*sigma./tmp(1:10);
%  t_knot1(8) = [];
% end
% 1/v(7,1)
% v/v(7,1)
% 
% return

S = 1;
[t_knot1, S, constraints] = FindStatisticallySignificantKnotRegions(t,x,Sigma,z_threshold,w, S);

[m_x1,Cm_x1,B1] = bspline_fit_no_tension_constrain(t,x,ones(size(x))*sigma,S,t_knot1,w, constraints);
tq1 = linspace(t(1),t(end),10*length(x))';
Bq1 = bspline(tq1,t_knot1,S+1);

tm1 = (t_knot1(1:end-1) + t_knot1(2:end))/2;
Bm1 = bspline(tm1,t_knot1,S+1);
val = squeeze(Bm1(:,:,2))*m_x1;
Ex = squeeze(Bm1(:,:,2))*Cm_x1*squeeze(Bm1(:,:,2)).';

x_fit1 = squeeze(Bq1(:,:,1))*m_x1;
x_fit1_error = sqrt(diag(squeeze(Bq1(:,:,1))*Cm_x1*squeeze(Bq1(:,:,1)).'));
x_error1 = x - squeeze(B1(:,:,1))*m_x1;
mean_x_error1 = sqrt(mean((path(tq1) - x_fit1).^2));
v_fit1 = squeeze(Bq1(:,:,2))*m_x1;
v_fit1_error = sqrt(diag(squeeze(Bq1(:,:,2))*Cm_x1*squeeze(Bq1(:,:,2)).'));
mean_v_error1 = sqrt(mean((speed(tq1) - v_fit1).^2));
mean_a_error1 = sqrt(mean((acceleration(tq0)).^2));

Q_error1 = (mean_x_error1 / sqrt(mean((path(tq1)).^2)) - 1) + (mean_v_error1 / sqrt(mean((speed(tq1)).^2)) - 1);

S = 2;
[t_knot2, S, constraints] = FindStatisticallySignificantKnotRegions(t,x,Sigma,z_threshold,w, S);
[m_x2,Cm_x2,B2] = bspline_fit_no_tension_constrain(t,x,ones(size(x))*sigma,S,t_knot2,w, constraints);
tq2 = linspace(t(1),t(end),10*length(x))';
Bq2 = bspline(tq2,t_knot2,S+1);


x_fit2 = squeeze(Bq2(:,:,1))*m_x2;
x_fit2_error = sqrt(diag(squeeze(Bq2(:,:,1))*Cm_x2*squeeze(Bq2(:,:,1)).'));
x_error2 = x - squeeze(B2(:,:,1))*m_x2;
mean_x_error2 = sqrt(mean((path(tq2) - x_fit2).^2));
v_fit2 = squeeze(Bq2(:,:,2))*m_x2;
v_fit2_error = sqrt(diag(squeeze(Bq2(:,:,2))*Cm_x2*squeeze(Bq2(:,:,2)).'));
mean_v_error2 = sqrt(mean((speed(tq2) - v_fit2).^2));
a_fit2 = squeeze(Bq2(:,:,3))*m_x2;
a_fit2_error = sqrt(diag(squeeze(Bq2(:,:,3))*Cm_x2*squeeze(Bq2(:,:,3)).'));
mean_a_error2 = sqrt(mean((acceleration(tq2) - a_fit2).^2));

Q_error2 = (mean_x_error2 / sqrt(mean((path(tq2)).^2)) - 1) + (mean_v_error2 / sqrt(mean((speed(tq2)).^2)) - 1) +(mean_a_error2 / sqrt(mean((acceleration(tq2)).^2)) - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Position Figures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

subplot(3,1,1)
PlotErrorRegion( tq0, x_fit0, z_threshold*x_fit0_error), hold on
plot(t,x_true,'k','LineWidth',1)
plot(tq0,x_fit0,'b','LineWidth',1.5)
scatter(t,x,6^2,'filled')
vlines(t_knot0, 'g--')
title(sprintf('Q-error: %.1f, position rms error=%.2f meters, velocity rms error=%.2f m/s, accel rms error=%.2f m/s^2', Q_error0,mean_x_error0, mean_v_error0, mean_a_error0))

subplot(3,1,2)
PlotErrorRegion( tq1, x_fit1, z_threshold*x_fit1_error), hold on
plot(t,x_true,'k','LineWidth',1)
plot(tq1,x_fit1,'b','LineWidth',1.5)
scatter(t,x,6^2,'filled')
vlines(t_knot1, 'g--')
title(sprintf('Q-error: %.1f, position rms error=%.2f meters, velocity rms error=%.2f m/s, accel rms error=%.2f m/s^2', Q_error1,mean_x_error1, mean_v_error1, mean_a_error1))

subplot(3,1,3)
PlotErrorRegion( tq2, x_fit2, z_threshold*x_fit2_error), hold on
plot(t,x_true,'k','LineWidth',1)
plot(tq2,x_fit2,'b','LineWidth',1.5)
scatter(t,x,6^2,'filled')
vlines(t_knot2, 'g--')
title(sprintf('Q-error: %.1f, position rms error=%.2f meters, velocity rms error=%.2f m/s, accel rms error=%.2f m/s^2', Q_error2,mean_x_error2, mean_v_error2, mean_a_error2))

figure
subplot(1,3,1)
hist(x_error0)
title(sprintf('std=%f',std(x_error0)))
subplot(1,3,2)
hist(x_error1)
title(sprintf('std=%f',std(x_error1)))
subplot(1,3,3)
hist(x_error2)
title(sprintf('std=%f',std(x_error2)))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Velocity Figures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[Diff1,t_v1,width] = FiniteDifferenceMatrixNoBoundary(1, t, 1);

figure

subplot(2,1,1)
PlotErrorRegion( tq1, v_fit1, z_threshold*v_fit1_error), hold on
plot(t,v_true,'k','LineWidth',1), hold on
plot( tq1, v_fit1,'b','LineWidth',1.5)
scatter(t_v1,Diff1*x,6^2,'filled')
vlines(t_knot1, 'g--')

subplot(2,1,2)
PlotErrorRegion( tq2, v_fit2, z_threshold*v_fit2_error), hold on
plot(t,v_true,'k','LineWidth',1)
plot( tq2, v_fit2,'b','LineWidth',1.5)
scatter(t_v1,Diff1*x,6^2,'filled')
vlines(t_knot2, 'g--')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Acceleration Figures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Diff2,t_a,width] = FiniteDifferenceMatrixNoBoundary(2, t, 1);

figure

PlotErrorRegion( tq2, a_fit2, z_threshold*a_fit2_error), hold on
plot(t,a_true,'k','LineWidth',1)
plot( tq2, a_fit2,'b','LineWidth',1.5)
scatter(t_a,Diff2*x)
vlines(t_knot2, 'g--')

% for i=2:length(t_knot)
%    range=(t_knot(i-1):t_knot(i))+1;
%    plot(t_a(range),mean(a(range))*ones(size(t_a(range))),'g')
% end
