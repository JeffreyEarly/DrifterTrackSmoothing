% Given some signal (t,x) contaminated by noise sigma, this uses the
% spectrum to estimate u_rms.
function [u_rms, u_std, u_mean] = EstimateRMSVelocityFromSpectrum( t, x, sigma, shouldPlotSpectra)

if length(unique(diff(t))) > 1
%    fprintf('interpolating...\n');
   dt = round(median(diff(t)));
   N = ceil((t(end)-t(1))/dt);
   t2 = dt*((0:(N-1))') + t(1);
   x = interp1(t,x,t2);
   t = t2;
end

[p,~,mu]=polyfit(t,x,1);
slope = p(1)/mu(2);
% intercept = p(2)-p(1)*mu(1)/mu(2);
u_mean = slope;

% now remove the linear trend
x = x-polyval(p,t,[],mu);

% first derivative, with points at t_u
[D,t_u] = FiniteDifferenceMatrixNoBoundary(1,t,1);
dt = t_u(2)-t_u(1);
T = t_u(end)-t_u(1);
N = length(t_u);

df = 1/T;
f = ([0:ceil(N/2)-1 -floor(N/2):-1]*df)';

ubar = fft(D*x);
s_signal = (dt/N)*ubar .* conj(ubar);

s_noise = sigma*sigma*dt*(2*pi*f).*(2*pi*f);

cutoff = 10;
% The factor of 10 is consitent with 80% confidence.
u2 = sum((s_signal > cutoff*s_noise) .* s_signal)*df;
u_std = sqrt(u2);
u_rms = sqrt( u2 + u_mean^2 );

if nargin > 3 && shouldPlotSpectra == 1
    f = fftshift(f);
    s_signal = fftshift(s_signal);
    s_noise = fftshift(s_noise);
    
    figure
    plot(f,s_signal)
    hold on
    plot(f,cutoff*s_noise), ylog
end