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
    if abs(meandiff) > 3*meansigma
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
    if abs(meandiff) > 3*meansigma
       a_knot_indices(end+1) = i; 
    end
end

t_knot1 = [t(1);  t(a_knot_indices+1); t(end)];

% for S=1:5
%     [Diff,~,width] = FiniteDifferenceMatrixNoBoundary(S, t, 1);
%     v = Diff*x;
%     Sigma=zeros(N-S,N-S);
%     for i=1:size(Sigma,1)
%         for j=1:size(Sigma,2)
%             Sigma(i,j) = sum(Diff(i,:).*Diff(j,:).*sigma'.*sigma');
%         end
%     end
%     Sigma_diag = diag(Sigma);
%     knot_indices = [];
%     %for i=2:(length(v)-S)
%     i=2;
%     while i<=length(v)-S
%         if isempty(knot_indices)
%             range = 1:(i-1);
%             mu = 0;
%         else
%             range=(knot_indices(end):(i-1));
%             mu = sum(width(range).*v(range))/sum(width(range));
%         end
%         meandiff = mean(v(i:i+S))-mu;
%         meansigma = sqrt(mean(Sigma_diag(range)));
%         if abs(meandiff) > 2.0*meansigma
%             knot_indices(end+1) = i;
%             i=i+1+S;
%         else
%             i=i+1;
%         end
%     end
%     
%     if isempty(knot_indices)
%        fprintf('Not significantly different from zero at order %d\n', S);
%        S=S-1;
%        break;
%     else
%         t_knot = [t(1);  t(knot_indices+1); t(end)];
%     end
% end

% for S=1:1
%     [Diff,~,width] = FiniteDifferenceMatrixNoBoundary(S, t, 1);
%     v = Diff*x;
%     Sigma=zeros(N-S,N-S);
%     for i=1:size(Sigma,1)
%         for j=1:size(Sigma,2)
%             Sigma(i,j) = sum(Diff(i,:).*Diff(j,:).*sigma'.*sigma');
%         end
%     end
%     Sigma_diag = diag(Sigma);
%     knot_indices = [];
%     for i=1:(length(v)-S)
%         if isempty(knot_indices)
%             range = i:(i+S);
%             sigmarange = range;
%             n_scale = 1/sqrt(S);
%             mu = 0;
%             meandiff = mean(v(range))-mu; % compare the mean velocity to zero.
%             meansigma = sqrt(sum(sum(Sigma(sigmarange,sigmarange),2),1));
%         else
%             range_prev = knot_indices(end):max((i-1),knot_indices(end)+S);
%             range_next = i:(i+S);
%             sigmarange = [range_prev,range_next];
%             n_scale = sqrt(1/length(range_prev) + 1/length(range_next));
%             meandiff = mean(v(range_next)) - mean(v(range_prev));
%             meansigma = sqrt(sum(sum(Sigma(range_prev,range_prev),2),1) + sum(sum(Sigma(range_next,range_next),2),1));
%         end
%         
%         
%         if abs(meandiff) > 2*meansigma*n_scale
%             knot_indices(end+1) = i;
%         end
%     end
%     
%     if isempty(knot_indices)
%        fprintf('Not significantly different from zero at order %d\n', S);
%         S=S-1;
%        break;
%     else
%         t_knot = [t(1);  t(knot_indices+1); t(end)];
%     end
% end
% S=S+1

for S=1:1
%     [Diff,~,width] = FiniteDifferenceMatrixNoBoundary(S, t, 1);
%     v = Diff*x;
%     Sigma=zeros(N-S,N-S);
%     for i=1:size(Sigma,1)
%         for j=1:size(Sigma,2)
%             Sigma(i,j) = sum(Diff(i,:).*Diff(j,:).*sigma'.*sigma');
%         end
%     end
   knot_indices = []; 
   v = [];
   Sigma_v = [];
   dt_v = [];
   
   for j=2:length(t)
        if isempty(knot_indices)
            i = 1;
            lastv = 0;
            lastsigma = 0;
            lastcov = 0;
        else
            i = knot_indices(end);
            lastv = v(end);
            lastsigma = Sigma_v(end);
            lastcov = sigma/dt_v(end);
        end
        
        newsigma = 2*sigma*sigma/( ( t(j) - t(i) )*( t(j) - t(i) ) );
        cov = -lastcov*sigma/( t(j) - t(i) );
        meansigma = sqrt( newsigma + lastsigma - 2*cov); % the sqrt(2) is the sample size in the z-test
        newv = ( x(j) - x(i) )/( t(j) - t(i) );
        if abs(newv-lastv) > 3*meansigma
            knot_indices(end+1) = j;
            v(end+1) = newv;
            Sigma_v(end+1) = newsigma;
            dt_v(end+1) = t(j) - t(i);
        end
   end
   
   t_knot = [t(1);  t(knot_indices); t(end)];
   
end

% Can we discern changes in position
knot_indices = [1];
t_knot = t(1);
Sigma = sigma*ones(size(t));
for j=2:length(t)
    i = knot_indices(end);
    range = i:j;
    
    dx = x(range)-mean(x(range));
    if (any(abs(dx) > 2.5*Sigma(range)))
        knot_indices(end+1) = j;
        t_knot(end+1) = t(j-1) + (t(j)-t(j-1))/2;
    end
end

if knot_indices(end) ~= length(t)
    knot_indices(end+1) = length(t);
    t_knot(end+1) = t(end);
end
t_knot = t_knot';

t_dx = t(knot_indices);
x_dx = x(knot_indices);
Sigma_dx = Sigma(knot_indices);
% At this point we'd have a piecewise constant spline. Each constant
% segment will average over points that are within 3\sigma of each other.
t_knot = t_dx;
S=1;

% knot_indices = [1];
% for j=2:length(t_dx)
%     i = knot_indices(end);
%     dx = x_dx(j)-x_dx(i);
%     s = sqrt( sigma(i)*sigma(i) + sigma(j)*sigma(j) );
%     if abs(dx) > 3*s
%         knot_indices(end+1) = j;
%     end
%     
%     if knot_indices(end) ~= length(t)
%         knot_indices(end+1) = length(t);
%     end
% end


Sigma = sigma*ones(size(t));
knot_indices = (1:length(t))';
left = knot_indices; % left most index of the grouping
right = knot_indices; % right most index of the grouping
xmean = zeros(size(left));  % mean of each grouping
for i=1:length(left)
    xmean(i) = mean(x(left(i):right(i)));
end
dx = diff(xmean); % difference between neighboring groupings

tolerance = zeros(size(dx));
for i=1:length(tolerance)
    tolerance(i) = sqrt(mean(Sigma(left(i):right(i+1)).^2));
end

% This isn't quite right. We need something based on a Chi-squared
% distribution. Because as our knowledge of the mean increase, we change
% our estimate of confidence.
z_score = abs(dx./tolerance);
[min_z_score,m_index] = min(z_score);

while (min_z_score < 2.5)
    % These two positions are indistinguishable, so merge them
    right(m_index) = right(m_index+1);
    left(m_index+1) = [];
    right(m_index+1) = [];
    
    xmean(m_index+1) = [];
    xmean(m_index) = mean(x(left(m_index):right(m_index)));
    
    dx(m_index) = [];
    tolerance(m_index) = [];
    if (m_index > 1) % update the left difference
        dx(m_index-1) = xmean(m_index)-xmean(m_index-1);
        tolerance(m_index-1) = sqrt(mean(Sigma(left(m_index-1):right(m_index)).^2));
    end
    if (m_index < length(xmean))
        dx(m_index) = xmean(m_index+1)-xmean(m_index);
        tolerance(m_index) = sqrt(mean(Sigma(left(m_index):right(m_index+1)).^2));
    end
    
    z_score = abs(dx./tolerance);
    [min_z_score,m_index] = min(z_score);
end

t_knot = [t(1);  (t(left(2:end))+t(right(1:(end-1))))/2; t(end)];

S=0;


[m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = drifter_fit_bspline_no_tension(t,x,x,ones(size(x))*sigma,ones(size(x))*sigma,S,t_knot,w);
x_fit = squeeze(Bq(:,:,1))*m_x;
% v_fit = squeeze(Bq(:,:,2))*m_x;
% a_fit = squeeze(Bq(:,:,3))*m_x;
% j_fit = squeeze(Bq(:,:,4))*m_x;
x_error = x - squeeze(B(:,:,1))*m_x;

% c(end+1)=a;
% d(end+1)=std(v_fit);
% f(end+1)=std(a_fit);
% g(end+1)=std(j_fit);

% end

mean_x_error = sqrt(mean((path(tq) - x_fit).^2));
% mean_v_error = sqrt(mean((speed(tq) - v_fit).^2));
mean_v_error=0;

figure
plot(t,x_true), hold on
plot(tq,x_fit,'b')
scatter(t,x,5)
title(sprintf('position rms error=%.2f meters, velocity rms error=%.2f m/s',mean_x_error, mean_v_error))

figure
hist(x_error)
title(sprintf('std=%f',std(x_error)))

return

figure
subplot(2,1,1)
plot( tq, v_fit,'b'), hold on
scatter(t_v,Diff1*x)
return
subplot(2,1,2)
plot( tq, a_fit,'b'), hold on
scatter(t_a,a)

% for i=2:length(t_knot)
%    range=(t_knot(i-1):t_knot(i))+1;
%    plot(t_a(range),mean(a(range))*ones(size(t_a(range))),'g')
% end
