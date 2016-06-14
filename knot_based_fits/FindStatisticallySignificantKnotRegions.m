function [t_knot, S, constraints] = FindStatisticallySignificantKnotRegions(t,x,Sigma,z_threshold,w,maxS)

knot_indices = (1:length(t))';
group = struct('left',knot_indices,'right',knot_indices,'value',[],'sigma2',[],'isNull',zeros(size(knot_indices)));
[t_knot, S, constraints] = GroupRecursion(0,group,t,x,Sigma,z_threshold, w, maxS);

end

% The names of this function are written as if group1 corresponds to the
% groupings for velocity, and group2 corresponds to the groupings for
% acceleration, but it is general and applies to all derivatives.
function group2 = CreateInitialGroupingForNextDerivative(S, group1)
group2 = struct('left',[],'right',[],'value',[],'sigma2',[],'isNull',[]);

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

group2.isNull = zeros(size(group2.left));

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

function constraints = ConstraintsFromGroup(group,t_knot,S)
constraints = struct('t',[],'K',[]);

for i=1:length(group.isNull)
    if group.isNull(i) == 1
        constraints.t(end+1) = (t_knot(i) + t_knot(i+1))/2;
        constraints.K(end+1) = S+1; 
    end
end

end

function [t_knot, S, constraints] = GroupRecursion(S,group,t,x,Sigma,z_threshold, w, maxS)

t_knot = KnotsFromGroup(group,t);
constraints = struct('t',[],'K',[]);
[merge_score, null_score, group] = ComputeZScore(S,group,t,x,Sigma,w,t_knot,constraints);
[min_z_score,m_index] = min(merge_score);
[min_null_score, m_null_index] = min(null_score);

while ( (~isempty(merge_score) && min_z_score < z_threshold) || (~isempty(null_score) && min_null_score < z_threshold))
    if (min_z_score < min_null_score) 
        % Merge groups that are not significantly different
        group.right(m_index) = group.right(m_index+1);
        group.isNull(m_index) = 0;
        group.left(m_index+1) = [];
        group.right(m_index+1) = [];
        group.isNull(m_index+1) = [];
    else
        group.isNull(m_null_index) = 1;
    end
    
    t_knot = KnotsFromGroup(group,t);
    constraints = ConstraintsFromGroup(group,t_knot,S);
    [merge_score, null_score, group] = ComputeZScore(S,group,t,x,Sigma,w,t_knot,constraints);
    [min_z_score,m_index] = min(merge_score);
    [min_null_score, m_null_index] = min(null_score);
end

if S == maxS
    return;
else
    group2 = CreateInitialGroupingForNextDerivative(S+1, group);
    [t_knot, S, constraints] = GroupRecursion(S+1,group2,t,x,Sigma,z_threshold, w, maxS);
    return;
end

end

function [merge_score, null_score, group] = ComputeZScore(S,group,t,x,Sigma,w,t_knot,constraints)

[m_x,Cm_x,~] = bspline_fit_no_tension_constrain(t,x,Sigma,S,t_knot,w,constraints);

t_error = (t_knot(1:end-1) + t_knot(2:end))/2;
B_error = bspline(t_error,t_knot,S+1);
group.value = squeeze(B_error(:,:,S+1))*m_x;
group.sigma2 = squeeze(B_error(:,:,S+1))*Cm_x*squeeze(B_error(:,:,S+1)).';

da = diff(group.value);
merge_tolerance = zeros(size(da));
for i=1:length(merge_tolerance)
    merge_tolerance(i) = sqrt(group.sigma2(i,i) + group.sigma2(i+1,i+1) - 2*group.sigma2(i,i+1));
end
merge_score = abs(da./merge_tolerance);

null_score = zeros(size(group.value));
for i=1:length(null_score)
    if group.isNull(i) == 1
        null_score(i) = Inf;
    else
        null_score(i) = abs(group.value(i)/sqrt(group.sigma2(i,i)));
    end
end

end