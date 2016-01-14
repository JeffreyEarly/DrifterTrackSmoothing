%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FindStatisticallySignificantChangesInVelocity
%
% This function returns knot at the point where there is a statistically
% significant change in velocity
%
% t             time vector, Nx1
% x             position vector, Nx1
% Sigma         error vector, Nx1
% z_threshold   cutoff for statistical significant (e.g., 3 => 3\sigma)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t_knot] = FindStatisticallySignificantChangesInVelocityFromGroup(group0,t,x,Sigma,z_threshold)

N = length(t);
S = 1;
[Diff,~,width] = FiniteDifferenceMatrixNoBoundary(S, t, 1);
v = Diff*x;
Sigma2x = Sigma.*Sigma;
Sigma2=zeros(N-S,N-S);
for i=1:size(Sigma2,1)
    for j=1:size(Sigma2,2)
        Sigma2(i,j) = sum(Diff(i,:).*Diff(j,:).*Sigma'.*Sigma');
    end
end

left = (1:(length(t)-1))'; % left most index of the grouping
right = (2:length(t))'; % right most index of the grouping
dv = diff(v); % difference between neighboring groupings

% the tolerance will be sqrt( sigma_left^2 + sigma_right^2 - 2*covariance)
tolerance = zeros(size(dv));
for i=1:length(tolerance)
    tolerance(i) = sqrt(Sigma2(i,i) + Sigma2(i+1,i+1) - 2*Sigma2(i,i+1));
end

z_score = abs(dv./tolerance);
[min_z_score,m_index] = min(z_score);

while (min_z_score < z_threshold)
    % These two positions are indistinguishable, so merge them
    right(m_index) = right(m_index+1);
    left(m_index+1) = [];
    right(m_index+1) = [];
    
    v(m_index+1) = [];
    v(m_index) =  (x(right(m_index)) - x(left(m_index)))/(t(right(m_index)) - t(left(m_index)));
    
    Sigma2(m_index+1,:) = [];
    Sigma2(:,m_index+1) = [];
    dv(m_index) = [];
    tolerance(m_index) = [];
    
    Sigma2(m_index,m_index) = (Sigma2x(left(m_index)) + Sigma2x(right(m_index)))/(t(right(m_index)) - t(left(m_index)))^2;
    if (m_index > 1) % update the left difference
        dv(m_index-1) = v(m_index)-v(m_index-1);
        cov = -Sigma2x(left(m_index))/( (t(right(m_index-1)) - t(left(m_index-1)))*(t(right(m_index)) - t(left(m_index))) );
        Sigma2(m_index-1,m_index) = cov;
        Sigma2(m_index,m_index-1) = cov;
        tolerance(m_index-1) = sqrt(Sigma2(m_index-1,m_index-1) + Sigma2(m_index,m_index) - 2*Sigma2(m_index-1,m_index));
    end
    if (m_index < length(v))
        dv(m_index) = v(m_index+1)-v(m_index);
        cov = -Sigma2x(left(m_index+1))/( (t(right(m_index)) - t(left(m_index)))*(t(right(m_index+1)) - t(left(m_index+1))) );
        Sigma2(m_index,m_index+1) = cov;
        Sigma2(m_index+1,m_index) = cov;
        tolerance(m_index) = sqrt(Sigma2(m_index,m_index) + Sigma2(m_index+1,m_index+1) - 2*Sigma2(m_index,m_index+1));
    end
    
    z_score = abs(dv./tolerance);
    [min_z_score,m_index] = min(z_score);
end

v_indices = [1; left(2:end); length(t)];
t_knot = t(v_indices);