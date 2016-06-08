function [t_knot, S, contraints] = FindStatisticallySignificantKnotRegions(t,x,Sigma,z_threshold,w)



end

function group2 = CreateGroupForNextDerivative(S, group1)
group2 = struct('left',[],'right',[],'value',[],'error',[]);

iAccelerationGroup = 1;
group2.left(iAccelerationGroup) = group1.left(1);

for iVelocityGroup = 1:length(group1.left)
    
    % if the velocity group ends spans more than one velocity, then...
    if (group1.right(iVelocityGroup) > group1.left(iVelocityGroup)+1 )
        group2.right(iAccelerationGroup) = group1.right(iVelocityGroup); % ...it's the end of this group
        
        if (iVelocityGroup == length(group1.left))
            continue;
        end
        
        iAccelerationGroup = iAccelerationGroup + 1;
        group2.left(iAccelerationGroup) = group1.right(iVelocityGroup)-1; % and the start of the next        
    end
    
    group2.right(iAccelerationGroup) = group1.right(iVelocityGroup)+1; % end the previous group
    
    iAccelerationGroup = iAccelerationGroup + 1;
    group2.left(iAccelerationGroup) = group1.right(iVelocityGroup); % and the start of the next
end



end

function [t_knot, S, contraints] = GroupRecursion(S,group,t,x,Sigma,z_threshold, w)



end