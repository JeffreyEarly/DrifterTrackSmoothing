% Return the difference between the max and the target acceleration
function totalError = AccelerationErrorForDespiking( sigma, nu, a, T, S, target_acceleration, drifters, shouldDisplay)

a = 10^(a);
w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));
Ndrifters = length(drifters.x);

accelerations = [];
for iDrifter = 1:Ndrifters
    x = drifters.x{iDrifter};
    y = drifters.y{iDrifter};
    t = drifters.t{iDrifter};
    
    tension = zeros(S,1);
    tension(T) = 1/a^2;
    [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = bspline_bivariate_fit_with_tension(t,x,y,ones(size(x))*sigma,ones(size(x))*sigma,S,tension, w);
    
    A = squeeze(B(:,:,3));
    ax = A*m_x;
    ay = A*m_y;
    
    accelerations = [accelerations; ax; ay];
end

totalError = abs(log10(target_acceleration/max(accelerations)));

if shouldDisplay == 1
    fprintf('(acceleration, lambda) = (%g, %f)\n', a, totalError);
end