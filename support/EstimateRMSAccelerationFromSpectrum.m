% Given some signal (t,x) contaminated by noise sigma, this uses the
% spectrum to estimate u_rms.
function [u_std, a_mean] = EstimateRMSAccelerationFromSpectrum( t, x, sigma)

if length(unique(diff(t))) > 1
   fprintf('interpolating...\n');
   dt = round(median(diff(t)));
   N = ceil((t(end)-t(1))/dt);
   t2 = dt*((0:(N-1))') + t(1);
   x = interp1(t,x,t2);
   t = t2;
end

[p,~,mu]=polyfit(t,x,2);
% slope = p(1)/mu(2);
% intercept = p(2)-p(1)*mu(1)/mu(2);
a_mean = 2*p(1)/mu(2)^2;

% now remove the quadratic trend
x = x-polyval(p,t,[],mu);

% first derivative, with points at t_u
[D,t_u] = FiniteDifferenceMatrixNoBoundary(2,t,1);
dt = t_u(2)-t_u(1);
T = t_u(end)-t_u(1);
nT = length(t_u);

fourierFrequencyT = 1/T;
f = ([0:ceil(nT/2)-1 -floor(nT/2):-1]*fourierFrequencyT)';

ubar = fft(D*x)*dt/sqrt(T);
s_signal = ubar .* conj(ubar);

s_noise = sigma*sigma*dt*(2*pi*f).^4;

u2 = sum((s_signal > 10.0*s_noise) .* s_signal)*fourierFrequencyT;
u_std = sqrt(u2);