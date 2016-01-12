%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FindStatisticallySignificantChangesInPosition
%
% This function returns knot at the space between statistically significant
% changes in position.
%
% t             time vector, Nx1
% x             position vector, Nx1
% Sigma         error vector, Nx1
% z_threshold   cutoff for statistical significant (e.g., 3 => 3\sigma)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t_knot] = FindStatisticallySignificantChangesInPosition(t,x,Sigma,z_threshold)

knot_indices = (1:length(t))';
left = knot_indices; % left most index of the grouping
right = knot_indices; % right most index of the grouping
xmean = zeros(size(left));  % mean of each grouping
for i=1:length(left)
    xmean(i) = mean(x(left(i):right(i)));
end
dx = diff(xmean); % difference between neighboring groupings

tolerance = zeros(size(dx));
for i=1:length(tolerance)
    tolerance(i) = sqrt(mean(Sigma(left(i):right(i+1)).^2));
end

% This isn't quite right. We need something based on a Chi-squared
% distribution. Because as our knowledge of the mean increase, we change
% our estimate of confidence.
z_score = abs(dx./tolerance);
[min_z_score,m_index] = min(z_score);

while (min_z_score < z_threshold)
    % These two positions are indistinguishable, so merge them
    right(m_index) = right(m_index+1);
    left(m_index+1) = [];
    right(m_index+1) = [];
    
    xmean(m_index+1) = [];
    xmean(m_index) = mean(x(left(m_index):right(m_index)));
    
    dx(m_index) = [];
    tolerance(m_index) = [];
    if (m_index > 1) % update the left difference
        dx(m_index-1) = xmean(m_index)-xmean(m_index-1);
        tolerance(m_index-1) = sqrt(mean(Sigma(left(m_index-1):right(m_index)).^2));
    end
    if (m_index < length(xmean))
        dx(m_index) = xmean(m_index+1)-xmean(m_index);
        tolerance(m_index) = sqrt(mean(Sigma(left(m_index):right(m_index+1)).^2));
    end
    
    z_score = abs(dx./tolerance);
    [min_z_score,m_index] = min(z_score);
end

t_knot = [t(1);  (t(left(2:end))+t(right(1:(end-1))))/2; t(end)];

