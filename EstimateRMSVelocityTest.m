addpath('./support');
load('sample_data/SyntheticTrajectories.mat')
sigma = position_error;

stride = 1;

if (stride < 10)
    shortenFactor = stride/10;
    DF = ceil(10/stride);
else
    shortenFactor = 1;
    DF = 1;
end

indices = 1:stride:floor(shortenFactor*length(t));
fprintf('Using %d points with stride %d\n', length(indices), stride);
x_obs = x(indices) + epsilon_x(indices);
y_obs = y(indices) + epsilon_y(indices);
t_obs = t(indices);

[D,t_u] = FiniteDifferenceMatrixNoBoundary(1,t_obs,1);
dt = t_u(2)-t_u(1);
T = t_u(end)-t_u(1);
N = length(t_u);

df = 1/T;
f = ([0:ceil(N/2)-1 -floor(N/2):-1]*df)';

u = D*x_obs;
v = D*y_obs;

ubar = fft(u)/N;
S_u_signal = (N*dt)* ubar .* conj(ubar);
vbar = fft(v)/N;
S_v_signal = (N*dt)* vbar .* conj(vbar);

u_noise = D*epsilon_x(indices);
v_noise = D*epsilon_y(indices);
ubar = fft(u_noise)/N;
S_u_noise = (N*dt)* ubar .* conj(ubar);
ubar = fft(v_noise)/N;
S_v_noise = (N*dt)* ubar .* conj(ubar);

mean((v_noise).^2)
sum(S_v_noise)*df

s_noise_estimate = sigma*sigma*dt*(2*pi*f).*(2*pi*f);
finitesum = (pi^2/3)*(sigma*sigma/dt^2)*(1+3/N + 2/N^2)
sum(s_noise_estimate)*df

% s_noise = sigma*sigma*dt*ones(size(f));


sum(s_noise_estimate)*df

std(D*epsilon_y(indices))^2
mean((D*epsilon_y(indices)).^2)
sum(s_noise_estimate)*df

cutoff = 15;
u_rms_estimate = sqrt(sum((S_u_signal > cutoff*s_noise_estimate) .* S_u_signal)*df);
v_rms_estimate = sqrt(sum((S_v_signal > cutoff*s_noise_estimate) .* S_v_signal)*df);
u_rms_true = sqrt(mean((D*x(indices)).^2));
v_rms_true = sqrt(mean((D*y(indices)).^2));

u_estimate_spectral = EstimateRMSVelocityFromSpectrum(t_obs,x_obs,sigma);
v_estimate_spectral = EstimateRMSVelocityFromSpectrum(t_obs,y_obs,sigma);

fprintf('u_rms_true: %f, u_rms_estimate: %f\n', u_rms_true, u_rms_estimate);
fprintf('v_rms_true: %f, v_rms_estimate: %f\n', v_rms_true, v_rms_estimate);

f = fftshift(f);
S_u_signal = fftshift(S_u_signal);
S_v_signal = fftshift(S_v_signal);
s_noise_estimate = fftshift(s_noise_estimate);

figure
plot(f,[S_u_signal, S_v_signal])
hold on
plot(f,s_noise_estimate), ylog