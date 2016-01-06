function [t_knot, S] = FindKnots(t, x, sigma, maxS)

N = length(t);

for S=1:maxS
    [Diff,~,width] = FiniteDifferenceMatrixNoBoundary(S, t, 1);
    v = Diff*x;
    Sigma=zeros(N-S,N-S);
    for i=1:size(Sigma,1)
        for j=1:size(Sigma,2)
            Sigma(i,j) = sum(Diff(i,:).*Diff(j,:).*sigma'.*sigma');
        end
    end
    Sigma_diag = diag(Sigma);
    knot_indices = [];
    for i=1:length(v)
        if isempty(knot_indices)
            range = 1:i;
            mu = 0;
        else
            range=(knot_indices(end):i);
            mu = sum(width(range).*v(range))/sum(width(range));
        end
        meandiff = v(i)-mu;
        meansigma = sqrt(mean(Sigma_diag(range)));
        if abs(meandiff) > 3*meansigma
            knot_indices(end+1) = i;
        end
    end
    
    if isempty(knot_indices)
       fprintf('Not significantly different from zero at order %d\n', S);
       S=S-1;
       return;
    else
        t_knot = [t(1);  t(knot_indices+1); t(end)];
    end
end


end