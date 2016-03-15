addpath('./support');
% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

load('sample_data/SyntheticTrajectories.mat')


result_stride = 2.^(0:8)';
result_rms_error = zeros(size(result_stride));
result_chi2_est = zeros(size(result_stride));
result_u_variance = zeros(size(result_stride));
result_a = zeros(size(result_stride));
result_dt = zeros(size(result_stride));
mean_standard_error = zeros(size(result_stride));
for i=1:length(result_stride)
    stride = result_stride(i);
    
    % Reduce the total length in some cases
    if (stride < 10)
        shortenFactor = stride/10;
        DF = ceil(10/stride);
    else
        shortenFactor = 1;
        DF = 1;
    end
    
    indices = 1:stride:floor(shortenFactor*length(t));
    fprintf('Using %d points with stride %d\n', length(indices), stride);
    indicesAll = 1:max(indices);
    x_obs = x(indices) + epsilon_x(indices);
    y_obs = y(indices) + epsilon_y(indices);
    t_obs = t(indices);
    sigma = position_error;

    S = 3;
    T = 2;
    K = S+1;
    a_vals = 10.^(linspace(-4,-2,15))';
    rms_error_T2 = zeros(size(a_vals));
    chi2_est_T2 = zeros(size(a_vals));
    final_variance_T2 = zeros(size(a_vals));
    standard_error_T2 = zeros(size(a_vals));
    for j = 1:length(a_vals)
       [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t_obs,x_obs,y_obs,sigma,S,T,a_vals(j), DF);
       X = squeeze(B(:,:,1));
       t_knot = NaturalKnotsForSpline( t_obs, K, DF );
       B = bspline(t(indicesAll),t_knot,K);
       X1 = squeeze(B(:,:,1));
       U1 = squeeze(B(:,:,T+1));
       rms_error_T2(j) = std( (x(indicesAll)-X1*m_x) + sqrt(-1)*(y(indicesAll)-X1*m_y) );
       chi2_est_T2(j) = std( (x_obs-X*m_x) + sqrt(-1)*(y_obs-X*m_y) );
       final_variance_T2(j) = std( U1*m_x + sqrt(-1)*U1*m_y )/sqrt(2);
       standard_error_T2(j) = (mean((diag(Cm_x))) + mean((diag(Cm_y))))/2;
    end
    
    [val, index] = min(rms_error_T2);
    fprintf('S=%d, T=2, stride=%d, min-a=%g, rms_error=%g, chi2-est=%g, a-variance=%g\n', S, stride, a_vals(index), val, (chi2_est_T2(index)/(position_error*sqrt(2)))^2, final_variance_T2(index) );
    result_rms_error(i) = rms_error_T2(index);
    result_chi2_est(i) = (chi2_est_T2(index)/(position_error*sqrt(2)))^2;
    result_u_variance(i) = final_variance_T2(index);
    result_a(i) = a_vals(index);
    result_dt(i) = t(result_stride(i)+1) - t(1);
    mean_standard_error(i) = standard_error_T2(index);
end
gamma = position_error./(0.20*result_dt);

figure
plot(1./gamma, result_chi2_est)
ylog

[p,S,mu]=polyfit(1./gamma,log10(result_chi2_est),1);
slope = p(1)/mu(2);
intercept = p(2)-p(1)*mu(1)/mu(2);

% This is the theoretical number of degrees of freedom that is being used.
nDOF = ((position_error*position_error)./mean_standard_error).*result_chi2_est

figure, plot( ((position_error*position_error)./mean_standard_error).*result_chi2_est, gamma)