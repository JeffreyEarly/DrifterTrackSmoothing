function [t_knot,group2] = FindStatisticallySignificantChangesInAccelFromGroupUsingRecu(group1,t,x,Sigma,z_threshold,w)

% we will use this as an index into group0
group2 = struct('left',[],'right',[],'value',[],'error',[]);

group2.left = []; % left most index of the grouping
group2.right = []; % right most index of the grouping

iAccelerationGroup = 1;
group2.left(iAccelerationGroup) = group1.left(1);
t_knot2 = t(1);
for iVelocityGroup = 1:length(group1.left)
    
    % if the velocity group ends spans more than one velocity, then...
    if (group1.right(iVelocityGroup) > group1.left(iVelocityGroup)+1 )
        group2.right(iAccelerationGroup) = group1.right(iVelocityGroup); % ...it's the end of this group
        
        if (iVelocityGroup == length(group1.left))
            continue;
        end
        
        iAccelerationGroup = iAccelerationGroup + 1;
        group2.left(iAccelerationGroup) = group1.right(iVelocityGroup)-1; % and the start of the next
        
        t_knot2(end+1) = ( t(group2.right(iAccelerationGroup-1)) + t(group2.left(iAccelerationGroup)) )/2;
    end
    
    group2.right(iAccelerationGroup) = group1.right(iVelocityGroup)+1; % end the previous group
    
    iAccelerationGroup = iAccelerationGroup + 1;
    group2.left(iAccelerationGroup) = group1.right(iVelocityGroup); % and the start of the next
    
    t_knot2(end+1) = ( t(group2.right(iAccelerationGroup-1)) + t(group2.left(iAccelerationGroup)) )/2;
    

end

group2.right(iAccelerationGroup) = group1.right(iVelocityGroup);
t_knot2(end+1) = t(end);

group2.left = group2.left(1:iAccelerationGroup)';
group2.right = group2.right(1:iAccelerationGroup)';
t_knot2 = t_knot2';

S = 2;
[m_x,Cm_x,B] = bspline_fit_no_tension(t,x,Sigma,S,t_knot2,w);

t_error = (t_knot2(1:end-1) + t_knot2(2:end))/2;
B_error = bspline(t_error,t_knot2,S+1);
group2.value = squeeze(B_error(:,:,S+1))*m_x;
group2.sigma2 = squeeze(B_error(:,:,S+1))*Cm_x*squeeze(B_error(:,:,S+1)).';

da = diff(group2.value);
tolerance = zeros(size(da));
for i=1:length(tolerance)
    tolerance(i) = sqrt(group2.sigma2(i,i) + group2.sigma2(i+1,i+1) - 2*group2.sigma2(i,i+1));
end

z_score = abs(da./tolerance);
[min_z_score,m_index] = min(z_score);

while (min_z_score < z_threshold)
    group2.right(m_index) = group2.right(m_index+1);
    group2.left(m_index+1) = [];
    group2.right(m_index+1) = [];
    
    t_knot = FindKnotsFromVelocityGroup(group2,t);
    
    [m_x,Cm_x,B] = bspline_fit_no_tension(t,x,Sigma,S,t_knot,w);
    
    t_error = (t_knot(1:end-1) + t_knot(2:end))/2;
    B_error = bspline(t_error,t_knot,S+1);
    group2.value = squeeze(B_error(:,:,S+1))*m_x;
    group2.sigma2 = squeeze(B_error(:,:,S+1))*Cm_x*squeeze(B_error(:,:,S+1))';
    
    da = diff(group2.value);
    tolerance = zeros(size(da));
    for i=1:length(tolerance)
        tolerance(i) = sqrt(group2.sigma2(i,i) + group2.sigma2(i+1,i+1) - 2*group2.sigma2(i,i+1));
    end
    
    z_score = abs(da./tolerance);
    [min_z_score,m_index] = min(z_score);
end

t_knot = FindKnotsFromVelocityGroup(group2,t);

% z_score = group2.value./sqrt(diag(group2.sigma2));
max_z_score = max(z_score);
fprintf('S2 Max z-score is %f\n', max_z_score);


end

function [t_knot] = FindKnotsFromVelocityGroup(group2,t)
t_knot = zeros(length(group2.left)+1,1);
t_knot(1) = t(1);
    for i=2:length(group2.left)
        t_knot(i)= (t(group2.left(i)) + t(group2.right(i-1)))/2;
    end
    t_knot(end) = t(end);
end