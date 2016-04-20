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
nT = length(t_u);

fourierFrequencyT = 1/T;
f = ([0:ceil(nT/2)-1 -floor(nT/2):-1]*fourierFrequencyT)';
df = f(2)-f(1);

ubar = fft(D*y_obs)/nT;
s_signal = (nT*dt)* ubar .* conj(ubar);

s_noise = sigma*sigma*dt*(2*pi*f).*(2*pi*f);
% s_noise = sigma*sigma*dt*ones(size(f));


sum(s_noise)*df

u2 = sum((s_signal > 10.0*s_noise) .* s_signal)*fourierFrequencyT;
u_rms_estimate = sqrt(u2);
u_rms_true = sqrt(mean((D*y(indices)).^2));

fprintf('u_rms_true: %f, u_rms_estimate: %f\n', u_rms_true, u_rms_estimate);
return;

figure
plot(f,s_signal)
hold on
plot(f,s_noise), ylog