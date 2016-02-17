reps = 100;
sigma = 300;
nu = 3.8;
t_noise = (0:2*60*60:20*86400)';
dt_noise = t_noise(2)-t_noise(1);
u_noise = reshape(StudentTNoise( sigma, nu, reps*length(t_noise)), [length(t_noise) reps]);
v_noise = reshape(StudentTNoise( sigma, nu, reps*length(t_noise)), [length(t_noise) reps]);
cv_noise = (diff(u_noise,1,1) + sqrt(-1)*diff(v_noise,1,1))/dt_noise;
[psi,lambda]=sleptap(size(cv_noise,1),taper_bandwidth);
[omega_noise,spp_noise,snn_noise,spn_noise]=mspec(dt_noise,cv_noise,psi);
f_noise=omega_noise*86400/(2*pi);

figure
plot(f_noise,vmean(spp_noise,2), 'black', 'LineWidth', 4), ylog, hold on,
plot(-f_noise,vmean(snn_noise,2), 'black', 'LineWidth', 4)