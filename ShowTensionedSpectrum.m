addpath('./support');
% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

load('sample_data/SyntheticTrajectories.mat')



stride = 1;

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

% Compute the spectrum
D = FiniteDifferenceMatrixNoBoundary(1,t(indices),1);
%  D = eye(length(t(indices)));
dt = t_obs(2)-t_obs(1);
cv = D*(x(indices) + sqrt(-1)*y(indices));
cepsilon = D*(epsilon_x(indices) + sqrt(-1)*epsilon_y(indices));
[~,spp,snn,spn]=mspec(dt,cv,[]);
[f,spp_e,snn_e,spn_e]=mspec(dt,cepsilon,[]);

df = f(2)-f(1);
(1/2/pi)*(sum(spp) + sum(snn(2:end)))*df
std(cv).^2 + mean(cv).*conj(mean(cv))

S = 3;
T = 2;
K = S+1;

t_knot = NaturalKnotsForSpline( t_obs, K, DF );
B = bspline(t(indicesAll),t_knot,K);
X1 = squeeze(B(:,:,1));
U1 = squeeze(B(:,:,T+1));

a = 5e-4;
errorFunction = @(a) TotalRMSError(t_obs,x_obs,y_obs,sigma,S,T, a, DF, X1, x(indicesAll), y(indicesAll));
optimalAcceleration = fminsearch( errorFunction, log10(a), optimset('TolX', 0.001, 'TolFun', 0.001) );
a = 10^(optimalAcceleration(1));

result_stride = zeros(3,1);
result_rms_error = zeros(size(result_stride));
result_chi2_est = zeros(size(result_stride));
result_u_variance = zeros(size(result_stride));
result_a = zeros(size(result_stride));
result_dt = zeros(size(result_stride));
mean_standard_error = zeros(size(result_stride));
result_snn = zeros(3,length(snn));
result_spp = zeros(3,length(spp));
result_snn_e = zeros(3,length(snn_e));
result_spp_e = zeros(3,length(spp_e));
result_s = zeros(3,length(snn_e));
result_s_e = zeros(3,length(spp_e));
for i=1:3
    result_a(i) = (2^(i-1))*a/2;
    
    [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t_obs,x_obs,y_obs,sigma,S,T,result_a(i), DF);
    X = squeeze(B(:,:,1));
    U = squeeze(B(:,:,2));
    
    result_rms_error(i) = std( (x(indicesAll)-X1*m_x) + sqrt(-1)*(y(indicesAll)-X1*m_y) );
    result_chi2_est(i) = (std( (x_obs-X*m_x) + sqrt(-1)*(y_obs-X*m_y) )/(position_error*sqrt(2)))^2;
    result_u_variance(i) = std( U1*m_x + sqrt(-1)*U1*m_y )/sqrt(2);
    
    result_dt(i) = t(result_stride(i)+1) - t(1);
    mean_standard_error(i) = (mean((diag(X*Cm_x*X.'))) + mean((diag(X*Cm_y*X.'))))/2;
    
    cv_obs = D*(X*m_x + sqrt(-1)*X*m_y);
    cepsilon_obs = D*((x_obs-X*m_x) + sqrt(-1)*(y_obs-X*m_y));
    [psi,lambda]=sleptap(size(cv_obs,1));
    [~,spp_obs,snn_obs,spn_obs]=mspec(dt,cv_obs,[]);
    [f,spp_e_obs,snn_e_obs,spn_e_obs]=mspec(dt,cepsilon_obs,[]);
    
    result_snn(i,:) = snn_obs;
    result_spp(i,:) = spp_obs;
    result_s(i,:) = vmean([snn_obs, spp_obs],2);
    result_snn_e(i,:) = snn_e_obs;
    result_spp_e(i,:) = spp_e_obs;
    result_s_e(i,:) = vmean([snn_e_obs, spp_e_obs],2);
    
    fprintf('S=%d, T=2, stride=%d, min-a=%g, rms_error=%g, chi2-est=%g, a-variance=%g\n', S, stride, a, result_rms_error(i), result_chi2_est(i), result_u_variance(i) );
end


timescale = 60;
ylimit = [1e-4 4e2];

figure
plot(f*timescale,sqrt(cumtrapz(f,snn+spp)/2/pi/2), 'k', 'LineWidth', 2)
hold on
plot(f*timescale,sqrt(cumtrapz(f,snn_e+spp_e)/2/pi/2), 'r', 'LineWidth', 2)
for i=1:3
    if i == 2
        lw = 2;
    else
        lw = 1;
    end
    plot(f*timescale,sqrt(cumtrapz(f,result_s(i,:))/2/pi),'LineWidth',lw)
    plot(f*timescale,sqrt(cumtrapz(f,result_s_e(i,:))/2/pi),'LineWidth',lw)
end
xlog, ylog
xlim([min(f*timescale) max(f*timescale)])


figure
plot(f*timescale,vmean([snn, spp],2), 'k', 'LineWidth', 2)
hold on
plot(f*timescale,vmean([snn_e, spp_e],2), 'r', 'LineWidth', 2)
for i=1:3
    if i == 2
        lw = 2;
    else
        lw = 1;
    end
    plot(f*timescale,result_s(i,:),'LineWidth',lw)
    plot(f*timescale,result_s_e(i,:),'LineWidth',lw)
end
xlog, ylog
xlim([min(f*timescale) max(f*timescale)])
ylim(ylimit)

figure
plot(f*timescale,vmean([snn + snn_e, spp + spp_e],2), 'k', 'LineWidth', 2)
hold on

for i=1:3
    if i == 2
        lw = 2;
    else
        lw = 1;
    end
    plot(f*timescale,result_s(i,:) + result_s_e(i,:),'LineWidth',lw)
end
xlog, ylog
xlim([min(f*timescale) max(f*timescale)])
ylim(ylimit)


df = f(2)-f(1);
s_total = vmean([snn + snn_e, spp + spp_e],2);
sum(s_total)*df/(2*pi)
sum(result_s + result_s_e,2)*df/(2*pi)

df = f(2)-f(1);
s_total = vmean([snn + snn_e, spp + spp_e],2);
sum((snn_e + spp_e)/2)*df/(2*pi)
sum(result_s_e,2)*df/(2*pi)

a = cumtrapz(f,result_s_e(2,:)/(2*pi) );
b = cumtrapz(f,(snn_e + spp_e)/2/(2*pi) );
figure, plot( f*timescale, a )
hold on, plot( f*timescale, b )

figure, plot(f*timescale,a'./b)

hold on, plot(f*timescale,a'./b, 'm')