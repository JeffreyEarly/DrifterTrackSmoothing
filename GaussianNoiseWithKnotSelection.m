% close all, clear

% c = []; d=[]; f=[]; g=[];



% Create a simple signal
dt_v = 1;
t=(0:dt_v:100)'; % seconds
a_true = 0.1;
u_true = 0; % meters/second

path = @(t) a_true*(t.*t) + u_true.*t;
speed = @(t) 2*a_true.*t + u_true*ones(size(t));

x_true = path(t);

sigma = 4; % meters

% Create some Gaussian noise
epsilon = sigma*randn(size(x_true));

% Contaminate the signal
x = x_true + epsilon;

N = length(t);

% Differentiation matrix for velocity, and velocity grid t_v
[Diff1,t_v] = FiniteDifferenceMatrixNoBoundary(1, t, 1);
v = Diff1*x;
dt_v = diff(t);
Sigma_x=zeros(N-1,N-1);
for i=1:size(Sigma_x,1)
    for j=1:size(Sigma_x,2)
        Sigma_x(i,j) = sum(Diff1(i,:).*Diff1(j,:).*sigma'.*sigma');
    end
end
sigma_x_diag = diag(Sigma_x);

v_knot_indices = [1]';
for i=2:(length(v)-1)
    range=(v_knot_indices(end):i);
    mu = sum(dt_v(range).*v(range))/sum(dt_v(range));
    meandiff = v(i)-mu;
    meansigma = sqrt(mean(sigma_x_diag(range)));
    if abs(meandiff) > 2*meansigma
       v_knot_indices(end+1) = i; 
    end
end
v_knot_indices(1) = [];

t_knot = [t(1); t(v_knot_indices); t(end)];

p = @(z) exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi));
w = @(z)(sigma*sigma);

S=2;
[Diff2,t_a] = FiniteDifferenceMatrixNoBoundary(2, t, 1);
a = Diff2*x;
dt_a = dt_v(1:end-1) + dt_v(2:end);
Sigma_xx=zeros(N-2,N-2);
for i=1:size(Sigma_xx,1)
    for j=1:size(Sigma_xx,2)
        Sigma_xx(i,j) = sum(Diff2(i,:).*Diff2(j,:).*sigma'.*sigma');
    end
end
sigma_xx_diag = diag(Sigma_xx);
a_knot_indices = [];
for i=1:length(a)
    if isempty(a_knot_indices)
        range = 1:i;
        mu = 0;
    else
        range=(a_knot_indices(end):i);
        mu = sum(dt_a(range).*a(range))/sum(dt_a(range));
    end
    meandiff = a(i)-mu;
    meansigma = sqrt(mean(sigma_xx_diag(range)));
    if abs(meandiff) > 2*meansigma
       a_knot_indices(end+1) = i; 
    end
end

t_knot = [t(1);  t(a_knot_indices); t(end)];



[m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = drifter_fit_bspline_no_tension(t,x,x,ones(size(x))*sigma,ones(size(x))*sigma,S,t_knot,w);
x_fit = squeeze(Bq(:,:,1))*m_x;
v_fit = squeeze(Bq(:,:,2))*m_x;
a_fit = squeeze(Bq(:,:,3))*m_x;
% j_fit = squeeze(Bq(:,:,4))*m_x;
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
plot(tq,x_fit,'b')
scatter(t,x,5)
title(sprintf('position rms error=%.2f meters, velocity rms error=%.2f m/s',mean_x_error, mean_v_error))

figure
hist(x_error)
title(sprintf('std=%f',std(x_error)))

figure
subplot(2,1,1)
plot( tq, v_fit,'b'), hold on
scatter(t_v,v)
subplot(2,1,2)
plot( tq, a_fit,'b'), hold on
scatter(t_a,a)