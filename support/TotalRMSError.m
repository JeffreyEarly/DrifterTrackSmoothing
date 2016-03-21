function rmsError = TotalRMSError(t_obs,x_obs,y_obs,sigma,S,T, a, DF, X1,x_all, y_all)
a = 10^a;
[m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t_obs,x_obs,y_obs,sigma,S,T, a, DF);
rmsError = std( (x_all-X1*m_x) + sqrt(-1)*(y_all-X1*m_y) );

shouldDisplay = 1;
if shouldDisplay == 1
    X = squeeze(B(:,:,1));
    SE = (mean((diag(X*Cm_x*X.'))) + mean((diag(X*Cm_y*X.'))))/2;
    fprintf('\t(acceleration, rms-error, standard-error) = (%g, %f, %f)\n', a, rmsError,SE);
end

end