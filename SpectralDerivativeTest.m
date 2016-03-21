addpath('./support');
% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

load('sample_data/SyntheticTrajectories.mat')

stride = 2;

% Reduce the total length in some cases
if (stride < 10)
    shortenFactor = stride/10;
    DF = ceil(10/stride);
else
    shortenFactor = 1;
    DF = 1;
end

% Reduce the number of points that we're considering
indices = 1:stride:floor(shortenFactor*length(t));
fprintf('Using %d points with stride %d\n', length(indices), stride);
indicesAll = 1:max(indices);
x_obs = x(indices) + epsilon_x(indices);
y_obs = y(indices) + epsilon_y(indices);
t_obs = t(indices);
sigma = position_error;


dt = t_obs(2)-t_obs(1);
cx = (x(indices) + sqrt(-1)*y(indices));
cepsilon = (epsilon_x(indices) + sqrt(-1)*epsilon_y(indices));


% Noise figure
figure
[f,spp_e,snn_e,spn_e]=mspec(dt,cepsilon,[]);
plot(f, vmean([spp_e, snn_e],2), 'k')
hold on
plot(f,2*sigma*sigma*dt*ones(size(f)),'k') % multiply by 2 b/c we want the one-sided spectrum
df = f(2)-f(1);
sum((snn_e + spp_e)/2)*df/(2*pi)
sum(2*sigma*sigma*dt*ones(size(f)))*df/(2*pi)

plot(f,f.*f.*vmean([spp_e, snn_e],2), 'k')
D = FiniteDifferenceMatrixNoBoundary(1,t(indices),1);
[f,spp_e,snn_e,spn_e]=mspec(dt,D*cepsilon,[]);
plot(f, vmean([spp_e, snn_e],2), 'r')
plot(f,2*sigma*sigma*dt*f.*f,'r')
legend('noise', 'noise theory', 'f*f*noise','d/dt noise', 'd/dt noise theory')
xlog, ylog

% Signal figure
figure
[psi,lambda]=sleptap(size(cx,1),10);
[f,spp_s,snn_s,spn_s]=mspec(dt,cx,psi);
plot(f, vmean([spp_s, snn_s],2), 'b', 'LineWidth', 2)
hold on
plot(f,sigma*sigma*dt*ones(size(f)),'k')
plot(f,f.*f.*vmean([spp_s, snn_s],2), 'g')
D = FiniteDifferenceMatrixNoBoundary(1,t(indices),1);
[f,spp_s,snn_s,spn_s]=mspec(dt,D*cx,[]);
plot(f, vmean([spp_s, snn_s],2), 'g', 'LineWidth', 2)
plot(f, vmean([spp_s, snn_s]./(f.*f),2), 'b', 'LineWidth', 2)
plot(f,sigma*sigma*dt*f.*f,'r')
xlog, ylog

cx_obs = (x_obs + sqrt(-1)*y_obs);
cv_obs = D*cx_obs;
[f,spp,snn,spn]=mspec(dt,cv_obs,[]);

signal = vmean([spp_e, snn_e]./(f.*f),2);
noise = 2*sigma*sigma*dt*ones(size(f));
noise_scalar = 2*sigma*sigma*dt;

signal = vmean([spp_e, snn_e],2);
noise = 2*sigma*sigma*dt.*f.*f;

figure
plot(f,signal)
hold on
plot(f,noise)
xlog,ylog

% Very important!!! You need to sum the signal, the measured thing, not the
% inverse of the measured thing, because that's what's unbiased.
ratio = signal ./ noise;

figure
plot(f, ratio, 'b', 'LineWidth', 2)
ylog

figure
plot(flip(f), cumsum(flip(ratio,1)))


plot(1:length(f), cumsum(flip(ratio,1)))