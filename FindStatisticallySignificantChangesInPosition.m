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

function [t_knot, group] = FindStatisticallySignificantChangesInPosition(t,x,Sigma,z_threshold)

group = struct('left',[],'right',[],'value',[]);

knot_indices = (1:length(t))';
group.left = knot_indices; % left most index of the grouping
group.right = knot_indices; % right most index of the grouping
group.value = zeros(size(group.left));  % mean of each grouping
for i=1:length(group.left)
    group.value(i) = mean(x(group.left(i):group.right(i)));
end
dx = diff(group.value); % difference between neighboring groupings

tolerance = zeros(size(dx));
for i=1:length(tolerance)
    tolerance(i) = sqrt(mean(Sigma(group.left(i):group.right(i+1)).^2)/length(group.left(i):group.right(i+1)));
end

% This isn't quite group.right. We need something based on a Chi-squared
% distribution. Because as our knowledge of the mean increase, we change
% our estimate of confidence.
z_score = abs(dx./tolerance);
[min_z_score,m_index] = min(z_score);

while (min_z_score < z_threshold)
    % These two positions are indistinguishable, so merge them
    group.right(m_index) = group.right(m_index+1);
    group.left(m_index+1) = [];
    group.right(m_index+1) = [];
    
    group.value(m_index+1) = [];
    group.value(m_index) = mean(x(group.left(m_index):group.right(m_index)));
    
    dx(m_index) = [];
    tolerance(m_index) = [];
    if (m_index > 1) % update the group.left difference
        dx(m_index-1) = group.value(m_index)-group.value(m_index-1);
        a = mean(Sigma(group.left(m_index-1):group.right(m_index-1)).^2)/length(group.left(m_index-1):group.right(m_index-1));
        b = mean(Sigma(group.left(m_index):group.right(m_index)).^2)/length(group.left(m_index):group.right(m_index));
        tolerance(m_index-1) = sqrt(a+b);
    end
    if (m_index < length(group.value))
        dx(m_index) = group.value(m_index+1)-group.value(m_index);
        a = mean(Sigma(group.left(m_index):group.right(m_index)).^2)/length(group.left(m_index):group.right(m_index));
        b = mean(Sigma(group.left(m_index+1):group.right(m_index+1)).^2)/length(group.left(m_index+1):group.right(m_index+1));
        tolerance(m_index) = sqrt(a+b);
    end
    
    z_score = abs(dx./tolerance);
    [min_z_score,m_index] = min(z_score);
end

t_knot = [t(1);  (t(group.left(2:end))+t(group.right(1:(end-1))))/2; t(end)];

