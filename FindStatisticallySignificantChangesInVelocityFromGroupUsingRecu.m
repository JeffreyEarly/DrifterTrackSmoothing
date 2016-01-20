function [t_knot] = FindStatisticallySignificantChangesInVelocityFromGroupUsingRecu(group0,t,x,Sigma,z_threshold,w)

% we will use this as an index into group0
group1 = struct('left',[],'right',[],'value',[],'error',[]);

group1.left = []; % left most index of the grouping
group1.right = []; % right most index of the grouping

iVelocityGroup = 0;
t_knot1 = [];
for iPositionGroup = 1:length(group0.left)
    
    if (iVelocityGroup > 0)
        group1.right(iVelocityGroup) = group0.left(iPositionGroup); % end the previous group
    end
    iVelocityGroup = iVelocityGroup + 1;
    group1.left(iVelocityGroup) = group0.left(iPositionGroup); % and the start of the next
    
    t_knot1(end+1) = t(group0.left(iPositionGroup));
    
    % if the position group ends on a different point, then...
    if (group0.right(iPositionGroup) ~= group0.left(iPositionGroup) )
        group1.right(iVelocityGroup) = group0.right(iPositionGroup); % ...it's the end of this group
        iVelocityGroup = iVelocityGroup + 1;
        group1.left(iVelocityGroup) = group0.right(iPositionGroup); % and the start of the next
        
        t_knot1(end+1) = t(group0.right(iPositionGroup));
    end
end
% last one doesn't count.
iVelocityGroup = iVelocityGroup - 1;

group1.left = group1.left(1:iVelocityGroup)';
group1.right = group1.right(1:iVelocityGroup)';
t_knot1 = t_knot1';

S = 1;
[m_x,Cm_x,B] = bspline_fit_no_tension(t,x,Sigma,S,t_knot1,w);

t_error = (t(group1.left)+t(group1.right))/2;
B_error = bspline(t_error,t_knot1,S+1);
group1.value = squeeze(B_error(:,:,S+1))*m_x;
group1.sigma2 = squeeze(B_error(:,:,S+1))*Cm_x*squeeze(B_error(:,:,S+1)).';

dv = diff(group1.value);
tolerance = zeros(size(dv));
for i=1:length(tolerance)
    tolerance(i) = sqrt(group1.sigma2(i,i) + group1.sigma2(i+1,i+1) - 2*group1.sigma2(i,i+1));
end

z_score = abs(dv./tolerance);
[min_z_score,m_index] = min(z_score);

while (min_z_score < z_threshold)
    group1.right(m_index) = group1.right(m_index+1);
    group1.left(m_index+1) = [];
    group1.right(m_index+1) = [];
    
    t_knot = FindKnotsFromVelocityGroup(group1);
    
    [m_x,Cm_x,B] = bspline_fit_no_tension(t,x,Sigma,S,t_knot,w);
    
    t_error = (t(group1.left)+t(group1.right))/2;
    B_error = bspline(t_error,t_knot,S+1);
    group1.value = squeeze(B_error(:,:,S+1))*m_x;
    group1.sigma2 = squeeze(B_error(:,:,S+1))*Cm_x*squeeze(B_error(:,:,S+1)).';
    
    dv = diff(group1.value);
    tolerance = zeros(size(dv));
    for i=1:length(tolerance)
        tolerance(i) = sqrt(group1.sigma2(i,i) + group1.sigma2(i+1,i+1) - 2*group1.sigma2(i,i+1));
    end
    
    z_score = abs(dv./tolerance);
    [min_z_score,m_index] = min(z_score);
end

end

function [t_knot] = FindKnotsFromVelocityGroup(group1)
t_knot = zeros(length(group1.left)+1,1);
    for i=1:length(group1.left)
        t_knot(i)=group1.left(i);
    end
    t_knot(end) = group1.right(end);
end