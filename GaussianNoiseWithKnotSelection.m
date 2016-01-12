% close all, clear

% c = []; d=[]; f=[]; g=[];



% Create a simple signal
dt_v = 1;
t=(0:dt_v:100)'; % seconds
a_true = 0.1;
u_true = 0; % meters/second

path = @(t) a_true*(t.*t) + u_true.*t;
speed = @(t) 2*a_true.*t + u_true*ones(size(t));

u_true = 500;
u_width = 10;
path = @(t) u_true*sech((t-50)/u_width).^2;
speed = @(t) -2*(u_true/u_width).*tanh((t-50)/u_width).*sech((t-50)/u_width).^2;

x_true = path(t);
v_true = speed(t);

sigma = 4; % meters

% Create some Gaussian noise
epsilon = sigma*randn(size(x_true));

% Contaminate the signal
x = x_true + epsilon;

N = length(t);

p = @(z) exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi));
w = @(z)(sigma*sigma);

Sigma = sigma*ones(size(t));
S = 0;
[t_knot] = FindStatisticallySignificantChangesInPosition(t,x,Sigma,3.0);
S = 1;
[t_knot] = FindStatisticallySignificantChangesInVelocity(t,x,Sigma,3.0);


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

[m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = drifter_fit_bspline_no_tension(t,x,x,ones(size(x))*sigma,ones(size(x))*sigma,S,t_knot,w);

x_fit = squeeze(Bq(:,:,1))*m_x;
x_error = x - squeeze(B(:,:,1))*m_x;
mean_x_error = sqrt(mean((path(tq) - x_fit).^2));

if S>0
    v_fit = squeeze(Bq(:,:,2))*m_x;
    mean_v_error = sqrt(mean((speed(tq) - v_fit).^2));
else
    mean_v_error=0;
end

figure
plot(t,x_true,'k','LineWidth',1), hold on
plot(tq,x_fit,'b','LineWidth',1.5)
scatter(t,x,6^2,'filled')
vlines(t_knot, 'g--')
title(sprintf('position rms error=%.2f meters, velocity rms error=%.2f m/s',mean_x_error, mean_v_error))

figure
hist(x_error)
title(sprintf('std=%f',std(x_error)))

if S>0
    [Diff1,t_v1,width] = FiniteDifferenceMatrixNoBoundary(S, t, 1);
    
    figure
    subplot(2,1,1)
    plot(t,v_true,'k','LineWidth',1), hold on
    plot( tq, v_fit,'b','LineWidth',1.5)
    scatter(t_v1,Diff1*x,6^2,'filled')
    vlines(t_knot, 'g--')
end

if S>1
    subplot(2,1,2)
    a_fit = squeeze(Bq(:,:,3))*m_x;
    plot( tq, a_fit,'b'), hold on
    scatter(t_a,Diff2*x)
    vlines(t_knot, 'g--')
end

% for i=2:length(t_knot)
%    range=(t_knot(i-1):t_knot(i))+1;
%    plot(t_a(range),mean(a(range))*ones(size(t_a(range))),'g')
% end
