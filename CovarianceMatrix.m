addpath('support');

t=(0:10)';
N = length(t);

a = [];
for S=1:7

% 2nd order 2nd derivative matrix, with 
% S = 4;20
D2 = FiniteDifferenceMatrixNoBoundary(S,t,1);

Sigma=zeros(N-S,N-S);
for i=1:size(Sigma,1)
    for j=1:size(Sigma,1)
        Sigma(i,j) = sum(D2(i,:).*D2(j,:));
    end
end

a = [a;Sigma(1,1)];

end

% D1 = FiniteDifferenceMatrixNoBoundary(1,t,1);
% D1(1:8,1:9)*D1(1:9,1:10)*D1

return;

dx = ones(size(t));

% 2nd order 2nd derivative matrix, no boundary terms
D2 = FiniteDifferenceMatrixNoBoundary(2,t,2);
D2 = D2(2:(N-1),:);

Sigma=zeros(N-2,N-2);
for i=1:size(Sigma,1)
    for j=1:size(Sigma,2)
        Sigma(i,j) = sum(D2(i,:).*D2(j,:).*dx'.*dx');
    end
end