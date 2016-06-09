function [t_knot, S, constraints] = FindStatisticallySignificantKnotRegions(t,x,Sigma,z_threshold,w,maxS)

knot_indices = (1:length(t))';
group = struct('left',knot_indices,'right',knot_indices,'value',[],'sigma2',[]);
[t_knot, S, constraints] = GroupRecursion(0,group,t,x,Sigma,z_threshold, w, maxS);

end

% The names of this function are written as if group1 corresponds to the
% groupings for velocity, and group2 corresponds to the groupings for
% acceleration, but it is general and applies to all derivatives.
function group2 = CreateInitialGroupingForNextDerivative(S, group1)
group2 = struct('left',[],'right',[],'value',[],'sigma2',[]);

iAccelerationGroup = 0;
for iVelocityGroup = 1:length(group1.left)
    % if the velocity group ends spans more than one velocity, then...
    if (group1.right(iVelocityGroup) > group1.left(iVelocityGroup)+S-1 )
        iAccelerationGroup = iAccelerationGroup + 1;
        group2.left(iAccelerationGroup) = group1.left(iVelocityGroup);
        group2.right(iAccelerationGroup) = group1.right(iVelocityGroup); % ...we keep the group as is...
        
        if (iVelocityGroup == length(group1.left)) % bail out of the for-loop if we've reached the end
            continue;
        end
        
        iAccelerationGroup = iAccelerationGroup + 1; % ...and create a transition group.
        group2.left(iAccelerationGroup) = group1.left(iVelocityGroup+1)-1; 
        group2.right(iAccelerationGroup) = group1.right(iVelocityGroup)+1; % can reach end, but not exceed
    elseif group1.right(iVelocityGroup)+1 > group1.right(end)
        continue;
    else
        iAccelerationGroup = iAccelerationGroup + 1;
        group2.left(iAccelerationGroup) = group1.left(iVelocityGroup);
        group2.right(iAccelerationGroup) = group1.right(iVelocityGroup)+1; %...otherwise grow the group by one        
    end

end

end

% This should be general
function [t_knot] = KnotsFromGroup(group,t)
t_knot = zeros(length(group.left)+1,1);
t_knot(1) = t(1);
    for i=2:length(group.left)
        t_knot(i)= (t(group.left(i)) + t(group.right(i-1)))/2;
    end
    t_knot(end) = t(end);
end

function [t_knot, S, constraints] = GroupRecursion(S,group,t,x,Sigma,z_threshold, w, maxS)

constraints = struct('t',[],'K',[]);

z_score = ComputeZScore(S,group,t,x,Sigma, w,constraints);
[min_z_score,m_index] = min(z_score);

[null_score, t_val] = ComputeNullScore(S,group,t,x,Sigma, w,constraints);
[min_null_score, m_null_index] = min(null_score);

while (min_z_score < z_threshold || min_null_score < z_threshold)
    if (min_z_score < min_null_score) 
        % Merge groups that are not significantly different
        group.right(m_index) = group.right(m_index+1);
        group.left(m_index+1) = [];
        group.right(m_index+1) = [];
    else
        constraints.t(end+1) = t_val(m_null_index);
        constraints.K(end+1) = S+1;
    end
    
    z_score = ComputeZScore(S,group,t,x,Sigma, w,constraints);
    [min_z_score,m_index] = min(z_score);
    
    [null_score, t_val] = ComputeNullScore(S,group,t,x,Sigma, w,constraints);
    [min_null_score, m_null_index] = min(null_score);
end

t_knot = KnotsFromGroup(group,t);

if S == maxS
    return;
else
    group2 = CreateInitialGroupingForNextDerivative(S+1, group);
    [t_knot, S, constraints] = GroupRecursion(S+1,group2,t,x,Sigma,z_threshold, w, maxS);
    return;
end

end

function [null_score, t_error] = ComputeNullScore(S,group,t,x,Sigma, w,constraints)
t_knot = KnotsFromGroup(group,t);

[m_x,Cm_x,~] = bspline_fit_no_tension_constrain(t,x,Sigma,S,t_knot,w,constraints);

t_error = (t_knot(1:end-1) + t_knot(2:end))/2;
B_error = bspline(t_error,t_knot,S+1);
group.value = squeeze(B_error(:,:,S+1))*m_x;
group.sigma2 = squeeze(B_error(:,:,S+1))*Cm_x*squeeze(B_error(:,:,S+1)).';

tolerance = zeros(size(group.value));
for i=1:length(tolerance)
    tolerance(i) = sqrt(group.sigma2(i,i));
end

null_score = abs(group.value./tolerance);
end

function z_score = ComputeZScore(S,group,t,x,Sigma, w,constraints)
t_knot = KnotsFromGroup(group,t);

[m_x,Cm_x,~] = bspline_fit_no_tension_constrain(t,x,Sigma,S,t_knot,w,constraints);

t_error = (t_knot(1:end-1) + t_knot(2:end))/2;
B_error = bspline(t_error,t_knot,S+1);
group.value = squeeze(B_error(:,:,S+1))*m_x;
group.sigma2 = squeeze(B_error(:,:,S+1))*Cm_x*squeeze(B_error(:,:,S+1)).';

da = diff(group.value);
tolerance = zeros(size(da));
for i=1:length(tolerance)
    tolerance(i) = sqrt(group.sigma2(i,i) + group.sigma2(i+1,i+1) - 2*group.sigma2(i,i+1));
end

z_score = abs(da./tolerance);
end