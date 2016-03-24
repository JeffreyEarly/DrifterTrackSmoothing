addpath('./support');
% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

load('sample_data/SyntheticTrajectories.mat')



result_stride = 2.^(0:8)';
% result_stride = [1;2;4;6;8;10;12;14;16;32;48;64;80;96;112;128;144;160;176;192;208;224;240;256];
result_rms_error = zeros(size(result_stride));
result_chi2_est = zeros(size(result_stride));
result_chi2_u_est = zeros(size(result_stride));
result_chi2_a_est = zeros(size(result_stride));
result_u_variance = zeros(size(result_stride));
result_a = zeros(size(result_stride));
result_dt = zeros(size(result_stride));
result_n = zeros(size(result_stride));
mean_standard_error = zeros(size(result_stride));
results_ac = zeros( 31, length(result_stride) );
results_ac_p = zeros( 30, length(result_stride) );
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

    D = FiniteDifferenceMatrixNoBoundary(1,t_obs,1);
    D2 = FiniteDifferenceMatrixNoBoundary(2,t_obs,1);
    
    S = 3;
    T = 2;
    K = S+1;
    a_vals = 10.^(linspace(-4,-2,15))';
    rms_error_T2 = zeros(size(a_vals));
    chi2_est_T2 = zeros(size(a_vals));
    final_variance_T2 = zeros(size(a_vals));
    standard_error_T2 = zeros(size(a_vals));
    
    t_knot = NaturalKnotsForSpline( t_obs, K, DF );
    B = bspline(t(indicesAll),t_knot,K);
    X1 = squeeze(B(:,:,1));
    U1 = squeeze(B(:,:,T+1));
    
    a = 5e-4;
    errorFunction = @(a) TotalRMSError(t_obs,x_obs,y_obs,sigma,S,T, a, DF, X1, x(indicesAll), y(indicesAll));
    optimalAcceleration = fminsearch( errorFunction, log10(a), optimset('TolX', 0.001, 'TolFun', 0.001) );
    a = 10^(optimalAcceleration(1));
    
    [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t_obs,x_obs,y_obs,sigma,S,T,a, DF);
    X = squeeze(B(:,:,1));
    
    dt = t(result_stride(i)+1) - t(1);
    results_ac(:,i) = Autocorrelation(x_obs-X*m_x,30);
    results_ac_p(:,i) =LjungBoxTest(results_ac(:,i),length(t_obs));
    result_rms_error(i) = std( (x(indicesAll)-X1*m_x) + sqrt(-1)*(y(indicesAll)-X1*m_y) );
    result_chi2_est(i) = (std( (x_obs-X*m_x) + sqrt(-1)*(y_obs-X*m_y) )/(position_error*sqrt(2)))^2;
    result_chi2_u_est(i) = ( std( D*((x_obs-X*m_x) + sqrt(-1)*(y_obs-X*m_y)) ) / ((sqrt(2)*position_error/dt)*sqrt(2)) )^2;
    result_chi2_a_est(i) = ( std( D2*((x_obs-X*m_x) + sqrt(-1)*(y_obs-X*m_y)) ) / ((sqrt(6)*position_error/dt^2)*sqrt(2)) )^2;
    result_u_variance(i) = std( U1*m_x + sqrt(-1)*U1*m_y )/sqrt(2);
    result_n(i) = length(t_obs);
    result_a(i) = a;
    result_dt(i) = dt;
    mean_standard_error(i) = (mean((diag(X*Cm_x*X.'))) + mean((diag(X*Cm_y*X.'))))/2;
    
    fprintf('S=%d, T=2, stride=%d, min-a=%g, rms_error=%g, chi2-est=%g, a-variance=%g\n', S, stride, a, result_rms_error(i), result_chi2_est(i), result_u_variance(i) );
end
gamma = position_error./(0.20*result_dt);

save('SampleVarianceVsGammaData.mat','result_rms_error', 'result_chi2_est','result_u_variance','result_a','result_dt','mean_standard_error','gamma', 'result_n');

gammaIndices = 1:24; gammaIndices(23) = [];
[p,S,mu]=polyfit(1./gamma(gammaIndices),log(result_chi2_est(gammaIndices)),1);
decay = 1/(p(1)/mu(2));
A = exp((p(2)-p(1)*mu(1)/mu(2)));

[p,S,mu]=polyfit(1./gamma(gammaIndices),(result_chi2_est(gammaIndices)),1);
m = 1/(p(1)/mu(2));
b = (p(2)-p(1)*mu(1)/mu(2));

figure
plot(1./gamma(gammaIndices), result_chi2_est(gammaIndices))
hold on, plot(1./gamma(gammaIndices), exp(- 1./(5*gamma(gammaIndices)) ), 'r');
plot(1./gamma(gammaIndices), - 1./(12*gamma(gammaIndices)) + 1 , 'g');
ylog



% This is the theoretical number of degrees of freedom that is being used.
gammaIndices = 1:24; gammaIndices(23) = [];
nDOF = ((position_error*position_error)./mean_standard_error);

[p,S,mu]=polyfit(log(gamma(gammaIndices)),log(nDOF(gammaIndices)-1),1);
slope = p(1)/mu(2)
intercept = p(2)-p(1)*mu(1)/mu(2)

[p,S,mu]=polyfit(1+gamma(gammaIndices),nDOF(gammaIndices),1);
slope = p(1)/mu(2);
intercept = p(2)-p(1)*mu(1)/mu(2);



figure, plot(gamma(gammaIndices), nDOF(gammaIndices))
hold on, plot(gamma(gammaIndices), slope*(1+gamma(gammaIndices)) + intercept);

gammaIndices = 1:22;
[p,S,mu]=polyfit(gamma(gammaIndices),(result_a(gammaIndices)/a_true).^2,1);
slope = p(1)/mu(2)
intercept = p(2)-p(1)*mu(1)/mu(2)

figure, plot(gamma(gammaIndices),(result_a(gammaIndices)/a_true).^2), hold on
plot(gamma(gammaIndices), slope*gamma(gammaIndices) + intercept)

nDOF2 = 1./(1-result_chi2_est);
[p,S,mu]=polyfit(log(gamma(gammaIndices)),log(nDOF2(gammaIndices)-1) ,1);
slope = p(1)/mu(2)
intercept = p(2)-p(1)*mu(1)/mu(2)
figure, plot(gamma(gammaIndices), nDOF2(gammaIndices)), hold on
plot(gamma(gammaIndices), slope*(1+gamma(gammaIndices)) + intercept)


[p,S,mu]=polyfit(log(nDOF(gammaIndices)-1),log(nDOF2(gammaIndices)-1) ,1);
slope = p(1)/mu(2)
intercept = p(2)-p(1)*mu(1)/mu(2)
