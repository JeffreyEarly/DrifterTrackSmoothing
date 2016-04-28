function error = DegreesOfFreedomError(t_obs,x_obs,y_obs,sigma,S,T, a, DF, expectedDOF)

a = 10^a;
[m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t_obs,x_obs,y_obs,sigma,S,T, a, DF);
X = squeeze(B(:,:,1));
SE = (mean((diag(X*Cm_x*X.'))) + mean((diag(X*Cm_y*X.'))))/2;
dof = sigma*sigma/SE;

error = abs(dof-expectedDOF);

shouldDisplay = 0;
if shouldDisplay == 1
    X = squeeze(B(:,:,1));
    SE = (mean((diag(X*Cm_x*X.'))) + mean((diag(X*Cm_y*X.'))))/2;
    fprintf('\t(acceleration, dof, expected-dof) = (%g, %f, %f)\n', a, dof,expectedDOF);
end

end