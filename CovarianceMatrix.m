t=(0:10)';
N = length(t);

% 2nd order 2nd derivative matrix, with 
D2 = FiniteDifferenceMatrix(2,t,2,2,2);

Sigma=zeros(size(D2));
for i=1:size(Sigma,1)
    for j=1:size(Sigma,2)
        Sigma(i,j) = sum(D2(i,:).*D2(j,:));
    end
end

dx = ones(size(t));

% 2nd order 2nd derivative matrix, no boundary terms
D2 = FiniteDifferenceMatrix(2,t,2,2,2);
D2 = D2(2:(N-1),:);

Sigma=zeros(N-2,N-2);
for i=1:size(Sigma,1)
    for j=1:size(Sigma,2)
        Sigma(i,j) = sum(D2(i,:).*D2(j,:).*dx'.*dx');
    end
end