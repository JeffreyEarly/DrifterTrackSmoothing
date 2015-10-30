function [D,z] = FiniteDifferenceMatrixNoBoundary(numDerivs, x, order)
% Creates a finite difference matrix of aribtrary accuracy, on an arbitrary
% grid. It does not implement boundary conditions (check my other routine
% for that), because it seeks to make all rows linearly independent.
%
% numDerivs ? the number of derivatives
% x ? the grid
% z location where approximations are to be accurate,
% orderOfAccuracy ? minimum order of accuracy required
%
% Jeffrey J. Early, 2015

n = length(x);
m = n - numDerivs;
D = zeros(m,n);

% order != accurracy.
nPoints = (numDerivs+1) + 2*(order-1);

if mod(numDerivs,2) == 0
    half = numDerivs/2;
    z = x((1+half):(n-half));
else
    mids = x(1:(n-1)) + diff(x)/2;
    half = floor(numDerivs/2);
    z = mids((1+half):(end-half));
end

% do we want to find the n closest points?
for i=1:m
    
    range_left = find( x <= z(i), ceil(nPoints/2), 'last');
    range_right = find( x > z(i), nPoints - length(range_left), 'first');
    range = union(range_left,range_right);
    
    if length(range)<nPoints
        range_right = find( x >= z(i), ceil(nPoints/2), 'first');
        range_left = find( x < z(i), nPoints - length(range_right), 'last');
        range = union(range_left,range_right);
    end
    
    c = weights( z(i), x(range), numDerivs );
    D(i,range) = c(numDerivs+1,:);
end

end

function c = weights(z,x,m)
% Calculates FD weights. The parameters are:
%  z   location where approximations are to be accurate,
%  x   vector with x-coordinates for grid points,
%  m   highest derivative that we want to find weights for
%  c   array size m+1,lentgh(x) containing (as output) in 
%      successive rows the weights for derivatives 0,1,...,m.
%
% Taken from Bengt Fornberg
%
    n=length(x); c=zeros(m+1,n); c1=1; c4=x(1)-z; c(1,1)=1;
    for i=2:n
       mn=min(i,m+1); c2=1; c5=c4; c4=x(i)-z;
       for j=1:i-1
          c3=x(i)-x(j);  c2=c2*c3;
          if j==i-1 
             c(2:mn,i)=c1*((1:mn-1)'.*c(1:mn-1,i-1)-c5*c(2:mn,i-1))/c2;
             c(1,i)=-c1*c5*c(1,i-1)/c2;
          end
          c(2:mn,j)=(c4*c(2:mn,j)-(1:mn-1)'.*c(1:mn-1,j))/c3;
          c(1,j)=c4*c(1,j)/c3;
       end
       c1=c2;
    end

end