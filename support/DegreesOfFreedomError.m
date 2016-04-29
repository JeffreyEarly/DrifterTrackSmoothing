function error = DegreesOfFreedomError(t_obs,x_obs,y_obs,sigma,S,T, log10lambda, DF, expectedDOF)

lambda = 10^log10lambda;
[m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t_obs,x_obs,y_obs,sigma,S,T, lambda, lambda, DF);
X = squeeze(B(:,:,1));
SE = (mean((diag(X*Cm_x*X.'))) + mean((diag(X*Cm_y*X.'))))/2;
dof = sigma*sigma/SE;

error = abs(dof-expectedDOF);

shouldDisplay = 0;
if shouldDisplay == 1
    X = squeeze(B(:,:,1));
    SE = (mean((diag(X*Cm_x*X.'))) + mean((diag(X*Cm_y*X.'))))/2;
    fprintf('\t(lambda, dof, expected-dof) = (%g, %f, %f)\n', lambda, dof,expectedDOF);
end

end