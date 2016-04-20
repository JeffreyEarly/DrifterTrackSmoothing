addpath('./support');
% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

load('sample_data/SyntheticTrajectories.mat')

result_stride = 2.^(0:8)';
result_stride = 1;

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
    
    % These are used to compute the rms-error using all data points.
    t_knot = NaturalKnotsForSpline( t_obs, K, DF );
    B = bspline(t(indicesAll),t_knot,K);
    X1 = squeeze(B(:,:,1));
    U1 = squeeze(B(:,:,T+1));

    u_estimate_spectral(i) = sqrt((EstimateRMSVelocityFromSpectrum(t_obs,x_obs,position_error)^2 + EstimateRMSVelocityFromSpectrum(t_obs,y_obs,position_error)^2)/2);
    a_estimate_spectral(i) = sqrt((EstimateRMSAccelerationFromSpectrum(t_obs,x_obs,position_error)^2 + EstimateRMSAccelerationFromSpectrum(t_obs,y_obs,position_error)^2)/2);
    
    dt = t_obs(2)-t_obs(1);
    expectedDOF(i) = 1 + 3*sigma/(u_estimate_spectral(i)*dt);
    
    % First blind guess as to what the optimal parameter should be
    a_blind_initial(i) = a_estimate_spectral(i)*sqrt(expectedDOF(i)/(expectedDOF(i)-1));
    
    % Now see how we did.
    [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t_obs,x_obs,y_obs,sigma,S,T,a_blind_initial(i), DF);
    X = squeeze(B(:,:,1));
    rms_error_blind_initial(i) = std( (x(indicesAll)-X1*m_x) + sqrt(-1)*(y(indicesAll)-X1*m_y) )/sqrt(2);
    dof_out_blind_initial(i) = sigma*sigma/( (mean((diag(X*Cm_x*X.'))) + mean((diag(X*Cm_y*X.'))))/2 );
    
    % Now optimize the parameter to match the expected DOF---still blind!
    errorFunction = @(a) DegreesOfFreedomError(t_obs,x_obs,y_obs,sigma,S,T, a, DF, expectedDOF(i));
    optimalAcceleration = fminsearch( errorFunction, log10(a_blind_initial(i)), optimset('TolX', 0.001, 'TolFun', 0.001) );
    a_blind_optimal(i) = 10^(optimalAcceleration(1));
    
    % Again, assess see how we did.
    [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t_obs,x_obs,y_obs,sigma,S,T,a_blind_optimal(i), DF);
    X = squeeze(B(:,:,1));
    rms_error_blind_optimal(i) = std( (x(indicesAll)-X1*m_x) + sqrt(-1)*(y(indicesAll)-X1*m_y) )/sqrt(2);
    dof_out_blind_optimal(i) = sigma*sigma/( (mean((diag(X*Cm_x*X.'))) + mean((diag(X*Cm_y*X.'))))/2 );
    
    % Finally, optimize the parameter to find the true best value we could
    % have set.
    errorFunction = @(a) TotalRMSError(t_obs,x_obs,y_obs,sigma,S,T, a, DF, X1, x(indicesAll), y(indicesAll));
    optimalAcceleration = fminsearch( errorFunction, log10(1e-4), optimset('TolX', 0.001, 'TolFun', 0.001) );
    a_true_optimal(i) = 10^(optimalAcceleration(1));
    
    [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t_obs,x_obs,y_obs,sigma,S,T,a_true_optimal(i), DF);
    X = squeeze(B(:,:,1));
    rms_error_true_optimal(i) = std( (x(indicesAll)-X1*m_x) + sqrt(-1)*(y(indicesAll)-X1*m_y) )/sqrt(2);
    dof_out_true_optimal(i) = sigma*sigma/( (mean((diag(X*Cm_x*X.'))) + mean((diag(X*Cm_y*X.'))))/2 );
    
    fprintf('S=%d, T=2, stride=%d, rms_error=%g, rms_error_blind_initial=%g, rms_error_blind_optimal=%g,\n', S, stride, rms_error_true_optimal(i), rms_error_blind_initial(i), rms_error_blind_optimal(i) );
end