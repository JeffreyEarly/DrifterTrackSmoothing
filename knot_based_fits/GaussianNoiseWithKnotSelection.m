% close all, clear

% c = []; d=[]; f=[]; g=[];

addpath('../support')


% Create a simple signal
dt_v = 1;
t=(0:dt_v:100)'; % seconds
a_true = 0.1;
u_true = 0; % meters/second

path = @(t) a_true*(t.*t) + u_true.*t;
speed = @(t) 2*a_true.*t + u_true*ones(size(t));

u_true = 500;
u_width = 10;
t_0 = 50;
path = @(t) u_true*sech((t-t_0)/u_width).^2;
speed = @(t) -2*(u_true/u_width).*tanh((t-t_0)/u_width).*sech((t-t_0)/u_width).^2;
acceleration = @(t) (u_true/u_width^2)*(4 * (tanh((t-t_0)/u_width)).^2 .* (sech((t-t_0)/u_width)).^2 - 2*(sech((t-t_0)/u_width)).^4);

x_true = path(t);
v_true = speed(t);
a_true = acceleration(t);

sigma = 200; % meters

% Create some Gaussian noise
epsilon = sigma*randn(size(x_true));

% Contaminate the signal
x = x_true + epsilon;

N = length(t);

p = @(z) exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi));
w = @(z)(sigma*sigma);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Now let's take a look at the acceleration and those errors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% v_indices = [1; left(2:end); length(t)];
% % v_indices = (1:length(t))';
% t_knot = t(v_indices);
% [Diff,t_v,width] = FiniteDifferenceMatrixNoBoundary(1, t_knot, 1);
% v = Diff*x(v_indices);
% Sigma2v = Sigma2;
% % S = 1;
% % Sigma2v=zeros(N-S,N-S);
% % for i=1:size(Sigma2v,1)
% %     for j=1:size(Sigma2v,2)
% %         Sigma2v(i,j) = sum(Diff(i,:).*Diff(j,:).*sigma'.*sigma');
% %     end
% % end
% 
% S = 2;
% [Diff,~,width] = FiniteDifferenceMatrixNoBoundary(2, t_knot, 1);
% a = Diff*x(v_indices);
% N = length(t_knot);
% Sigma2x = sigma*sigma*ones(N,1);
% Sigma2=zeros(N-S,N-S);
% for i=1:size(Sigma2,1)
%     for j=1:size(Sigma2,2)
%         Sigma2(i,j) = sum(Diff(i,:).*Diff(j,:).*sigma'.*sigma');
%     end
% end
% 
% left = (1:(length(v_indices)-1))'; % left most index of the grouping
% right = (2:length(v_indices))'; % right most index of the grouping
% da = diff(a); % difference between neighboring groupings
% 
% % the tolerance will be sqrt( sigma_left^2 + sigma_right^2 - 2*covariance)
% tolerance = zeros(size(da));
% for i=1:length(tolerance)
%     tolerance(i) = sqrt(Sigma2(i,i) + Sigma2(i+1,i+1) - 2*Sigma2(i,i+1));
% end
% 
% z_score = abs(da./tolerance);
% [min_z_score,m_index] = min(z_score);
% 
% while (min_z_score < 3.0)
%     % These two accelerations are indistinguishable, so merge them
%     right(m_index) = right(m_index+1);
%     left(m_index+1) = [];
%     right(m_index+1) = [];
%     
%     a(m_index+1) = [];
%     a(m_index) =  (v(right(m_index)) - v(left(m_index)))/(t_v(right(m_index)) - t_v(left(m_index)));
%     
%     Sigma2(m_index+1,:) = [];
%     Sigma2(:,m_index+1) = [];
%     da(m_index) = [];
%     tolerance(m_index) = [];
%     
%     Sigma2(m_index,m_index) = (Sigma2v(left(m_index),left(m_index)) + Sigma2v(right(m_index),right(m_index)))/(t_v(right(m_index)) - t_v(left(m_index)))^2;
%     if (m_index > 1) % update the left difference
%         da(m_index-1) = a(m_index)-a(m_index-1);
%         cov = -Sigma2v(left(m_index),left(m_index))/( (t_v(right(m_index-1)) - t_v(left(m_index-1)))*(t_v(right(m_index)) - t_v(left(m_index))) );
%         Sigma2(m_index-1,m_index) = cov;
%         Sigma2(m_index,m_index-1) = cov;
%         tolerance(m_index-1) = sqrt(Sigma2(m_index-1,m_index-1) + Sigma2(m_index,m_index) - 2*Sigma2(m_index-1,m_index));
%     end
%     if (m_index < length(a))
%         da(m_index) = a(m_index+1)-a(m_index);
%         cov = -Sigma2v(left(m_index+1),left(m_index+1))/( (t_v(right(m_index)) - t_v(left(m_index)))*(t_v(right(m_index+1)) - t_v(left(m_index+1))) );
%         Sigma2(m_index,m_index+1) = cov;
%         Sigma2(m_index+1,m_index) = cov;
%         tolerance(m_index) = sqrt(Sigma2(m_index,m_index) + Sigma2(m_index+1,m_index+1) - 2*Sigma2(m_index,m_index+1));
%     end
%     
%     z_score = abs(da./tolerance);
%     [min_z_score,m_index] = min(z_score);
% end
% 
% t_knot = [t(1); mean( [t(v_indices(left(2:(end-1)))), t(v_indices(right(2:(end-1))))],2 ); t(end)];

z_threshold = 4;

Sigma = sigma*ones(size(t));
S = 0;
[t_knot0, group0] = FindStatisticallySignificantChangesInPosition(t,x,Sigma,z_threshold);
[m_x0,Cm_x0,B0] = bspline_fit_no_tension(t,x,Sigma,S,t_knot0,w);

tq0 = linspace(t(1),t(end),10*N)';
Bq0 = bspline(tq0,t_knot0,S+1);
tm0 = (t(group0.left)+t(group0.right))/2;
Bm0 = bspline(tm0,t_knot0,S+1);
val = squeeze(Bm0(:,:,1))*m_x0;
Ex = squeeze(Bm0(:,:,1))*Cm_x0*squeeze(Bm0(:,:,1)).';
% A == B
% A = diag(Ex);
% B = group0.error;

x_fit0 = squeeze(Bq0(:,:,1))*m_x0;
x_error0 = x - squeeze(B0(:,:,1))*m_x0;
mean_x_error0 = sqrt(mean((path(tq0) - x_fit0).^2));

[t_knot1, group1] = FindStatisticallySignificantChangesInVelocityFromGroupUsingRecu(group0,t,x,Sigma,z_threshold,w);

S = 1;

% [m_x1,m_y1,Cm_x1,Cm_y1,B1,Bq1,tq1] = drifter_fit_bspline_no_tension(t,x,x,ones(size(x))*sigma,ones(size(x))*sigma,S,t_knot1,w);
[m_x1,Cm_x1,B1] = bspline_fit_no_tension(t,x,ones(size(x))*sigma,S,t_knot1,w);
tq1 = linspace(t(1),t(end),10*length(x))';
Bq1 = bspline(tq1,t_knot1,S+1);

tm1 = (t_knot1(1:end-1) + t_knot1(2:end))/2;
Bm1 = bspline(tm1,t_knot1,S+1);
val = squeeze(Bm1(:,:,2))*m_x1;
Ex = squeeze(Bm1(:,:,2))*Cm_x1*squeeze(Bm1(:,:,2)).';

x_fit1 = squeeze(Bq1(:,:,1))*m_x1;
x_error1 = x - squeeze(B1(:,:,1))*m_x1;
mean_x_error1 = sqrt(mean((path(tq1) - x_fit1).^2));
v_fit1 = squeeze(Bq1(:,:,2))*m_x1;
mean_v_error1 = sqrt(mean((speed(tq1) - v_fit1).^2));

[t_knot2, group2] = FindStatisticallySignificantChangesInAccelFromGroupUsingRecu(group1,t,x,Sigma,z_threshold,w);

S = 2;
% [m_x2,m_y2,Cm_x2,Cm_y2,B2,Bq2,tq2] = drifter_fit_bspline_no_tension(t,x,x,ones(size(x))*sigma,ones(size(x))*sigma,S,t_knot2,w);
[m_x2,Cm_x2,B2] = bspline_fit_no_tension(t,x,ones(size(x))*sigma,S,t_knot2,w);
tq2 = linspace(t(1),t(end),10*length(x))';
Bq2 = bspline(tq2,t_knot2,S+1);

x_fit2 = squeeze(Bq2(:,:,1))*m_x2;
x_error2 = x - squeeze(B2(:,:,1))*m_x2;
mean_x_error2 = sqrt(mean((path(tq2) - x_fit2).^2));
v_fit2 = squeeze(Bq2(:,:,2))*m_x2;
mean_v_error2 = sqrt(mean((speed(tq2) - v_fit2).^2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Position Figures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

subplot(3,1,1)
plot(t,x_true,'k','LineWidth',1), hold on
plot(tq0,x_fit0,'b','LineWidth',1.5)
scatter(t,x,6^2,'filled')
vlines(t_knot0, 'g--')
title(sprintf('position rms error=%.2f meters',mean_x_error0))

subplot(3,1,2)
plot(t,x_true,'k','LineWidth',1), hold on
plot(tq1,x_fit1,'b','LineWidth',1.5)
scatter(t,x,6^2,'filled')
vlines(t_knot1, 'g--')
title(sprintf('position rms error=%.2f meters, velocity rms error=%.2f m/s',mean_x_error1, mean_v_error1))

subplot(3,1,3)
plot(t,x_true,'k','LineWidth',1), hold on
plot(tq2,x_fit2,'b','LineWidth',1.5)
scatter(t,x,6^2,'filled')
vlines(t_knot2, 'g--')
title(sprintf('position rms error=%.2f meters, velocity rms error=%.2f m/s',mean_x_error2, mean_v_error2))


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
plot(t,v_true,'k','LineWidth',1), hold on
plot( tq1, v_fit1,'b','LineWidth',1.5)
scatter(t_v1,Diff1*x,6^2,'filled')
vlines(t_knot1, 'g--')

subplot(2,1,2)
plot(t,v_true,'k','LineWidth',1), hold on
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

a_fit2 = squeeze(Bq2(:,:,3))*m_x2;
plot(t,a_true,'k','LineWidth',1), hold on,
plot( tq2, a_fit2,'b','LineWidth',1.5)
scatter(t_a,Diff2*x)
vlines(t_knot2, 'g--')

% for i=2:length(t_knot)
%    range=(t_knot(i-1):t_knot(i))+1;
%    plot(t_a(range),mean(a(range))*ones(size(t_a(range))),'g')
% end
