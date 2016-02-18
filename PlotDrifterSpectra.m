addpath('./support');


% Choose which drifters to analyze.
load('sample_data/smoothedGriddedRho1Drifters.mat');
% load('griddedRho1Drifters.mat');

% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 2;
LoadFigureDefaults

numDrifters = size(x,2);
days = t/86400;
dt = t(2)-t(1);

f0 = abs(corfreq(lat0)/3600);

[x_com, y_com, q, r] = CenterOfMass( x, y );
cv_orig = (diff(x,1,1) + sqrt(-1)*diff(y,1,1))/dt;
cv = (diff(q,1,1) + sqrt(-1)*diff(r,1,1))/dt;

averaging_bandwidth=4;
taper_bandwidth=4;


[psi,lambda]=sleptap(size(cv,1),taper_bandwidth);
[~,spp_orig,snn_orig,spn_orig]=mspec(dt,cv_orig,psi);
[omega_p,spp,snn,spn]=mspec(dt,cv,psi);


% convert from radians/second to cycles/day
f_p=omega_p*86400/(2*pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create some noise!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reps = 100;
sigma = 8.0;
nu = 5.5;
t_noise = t(1:2:end);
dt_noise = t_noise(2)-t_noise(1);
u_noise = reshape(StudentTNoise( sigma, nu, reps*length(t_noise)), [length(t_noise) reps]);
v_noise = reshape(StudentTNoise( sigma, nu, reps*length(t_noise)), [length(t_noise) reps]);
cv_noise = (diff(u_noise,1,1) + sqrt(-1)*diff(v_noise,1,1))/dt_noise;
[psi,lambda]=sleptap(size(cv_noise,1),taper_bandwidth);
[omega_noise,spp_noise,snn_noise,spn_noise]=mspec(dt_noise,cv_noise,psi);
f_noise=omega_noise*86400/(2*pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Units', 'points', 'Position', [50 50 figure_width_2col 300*scaleFactor])
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

% Plot the individual drifters
plot(f_p,spp, 'Color', [0.4 0.4 0.4], 'LineWidth', 2),ylog
hold on
plot(-f_p,snn, 'Color', [0.4 0.4 0.4], 'LineWidth', 2)

% Plot the mean before strain removal
plot(f_p,vmean(spp_orig,2), 'green', 'LineWidth', 4)
plot(-f_p,vmean(snn_orig,2), 'green', 'LineWidth', 4)

% Plot the mean after strain removal
plot(f_p,vmean(spp,2), 'blue', 'LineWidth', 4)
plot(-f_p,vmean(snn,2), 'blue', 'LineWidth', 4)


% Plot the mean after strain removal
plot(f_noise,vmean(spp_noise,2), 'black', 'LineWidth', 4)
plot(-f_noise,vmean(snn_noise,2), 'black', 'LineWidth', 4)

xlim([-24 24])
%     min_s = 5e-2;
% 	max_s = 3e1;
%  	ylim([min_s max_s])

set(gca,'FontSize',figure_axis_tick_size)
xlabel('frequency (cycles per day)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('power (m^2/s)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
