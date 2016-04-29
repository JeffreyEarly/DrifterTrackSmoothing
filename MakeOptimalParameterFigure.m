addpath('./support');
% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

load('sample_data/SyntheticTrajectories.mat')

% Do you want to assess the error using all the points from the signal
% (which makes sense for an interpolation based metric) or just points from
% the observed (decimated) signal only?
shouldUseObservedSignalOnly = 1;

result_stride = 2.^(0:9)';
% result_stride = 1;

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
    if (shouldUseObservedSignalOnly == 1)
        indicesAll = indices;
    else
        indicesAll = 1:max(indices);
    end
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
    lambda_blind_initial(i) = (expectedDOF(i)-1)/(expectedDOF(i)*a_estimate_spectral(i)^2);
    
    % Now see how we did.
    [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t_obs,x_obs,y_obs,sigma,S,T,lambda_blind_initial(i),lambda_blind_initial(i), DF);
    X = squeeze(B(:,:,1));
    rms_error_blind_initial(i) = std( (x(indicesAll)-X1*m_x) + sqrt(-1)*(y(indicesAll)-X1*m_y) )/sqrt(2);
    dof_out_blind_initial(i) = sigma*sigma/( (mean((diag(X*Cm_x*X.'))) + mean((diag(X*Cm_y*X.'))))/2 );
    
    % Now optimize the parameter to match the expected DOF---still blind!
    errorFunction = @(log10lambda) DegreesOfFreedomError(t_obs,x_obs,y_obs,sigma,S,T,log10lambda,DF,expectedDOF(i));
    optimalLog10lambda = fminsearch( errorFunction, log10(lambda_blind_initial(i)), optimset('TolX', 0.001, 'TolFun', 0.001) );
    lambda_blind_optimal(i) = 10^(optimalLog10lambda);
    
    % Again, assess see how we did.
    [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t_obs,x_obs,y_obs,sigma,S,T,lambda_blind_optimal(i),lambda_blind_optimal(i),DF);
    X = squeeze(B(:,:,1));
    rms_error_blind_optimal(i) = std( (x(indicesAll)-X1*m_x) + sqrt(-1)*(y(indicesAll)-X1*m_y) )/sqrt(2);
    dof_out_blind_optimal(i) = sigma*sigma/( (mean((diag(X*Cm_x*X.'))) + mean((diag(X*Cm_y*X.'))))/2 );
    
    % Finally, optimize the parameter to find the true best value we could
    % have set.
    errorFunction = @(log10lambda) TotalRMSError(t_obs,x_obs,y_obs,sigma,S,T, log10lambda, DF, X1, x(indicesAll), y(indicesAll));
    optimalLog10lambda = fminsearch( errorFunction, log10(lambda_blind_optimal(i)), optimset('TolX', 0.001, 'TolFun', 0.001) );
    lambda_true_optimal(i) = 10^(optimalLog10lambda);
    
    [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t_obs,x_obs,y_obs,sigma,S,T,lambda_true_optimal(i), lambda_true_optimal(i), DF);
    X = squeeze(B(:,:,1));
    rms_error_true_optimal(i) = std( (x(indicesAll)-X1*m_x) + sqrt(-1)*(y(indicesAll)-X1*m_y) )/sqrt(2);
    dof_out_true_optimal(i) = sigma*sigma/( (mean((diag(X*Cm_x*X.'))) + mean((diag(X*Cm_y*X.'))))/2 );
    
    fprintf('S=%d, T=2, stride=%d, rms_error=%g, rms_error_blind_initial=%g, rms_error_blind_optimal=%g,\n', S, stride, rms_error_true_optimal(i), rms_error_blind_initial(i), rms_error_blind_optimal(i) );
end

if (shouldUseObservedSignalOnly == 1)
    outputFile = 'OptimalParametersObservedOnly.mat';
else
    outputFile = 'OptimalParameters.mat';
end

save(outputFile, 'u_estimate_spectral', 'a_estimate_spectral', 'expectedDOF', 'lambda_blind_initial', 'rms_error_blind_initial', 'dof_out_blind_initial', 'lambda_blind_optimal', 'rms_error_blind_optimal', 'dof_out_blind_optimal', 'lambda_true_optimal', 'rms_error_true_optimal', 'dof_out_true_optimal', 'result_stride')

for i=1:length(result_stride)
   fprintf('%d & %#.3g m (%#.3g) &  %#.3g m (%#.3g) &  %#.3g m (%#.3g) \\\\ \n', result_stride(i), rms_error_true_optimal(i), dof_out_true_optimal(i), rms_error_blind_optimal(i), dof_out_blind_optimal(i), rms_error_blind_initial(i), dof_out_blind_initial(i) )  ;
end