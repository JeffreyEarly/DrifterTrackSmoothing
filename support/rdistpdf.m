function [qq] = rdistpdf(nu,delta,mx)
% numerically computed PDF of r=\sqrt{x^+y^2}, where x and y are
% independent t(nu) distributions
% delta: increments in the pdf vector
% mx: maximum value of the pdf
% output vector qq is the pdf at locations delta:delta:mx
x=delta:delta:mx;
pp=zeros(1,1/delta*mx); qq=zeros(1,1/delta*mx);
for ii = 1:1/delta*mx
    gridX = x(ii)^2/10000;
    y=gridX/2:gridX:x(ii)^2-gridX/2;
    pp(ii)=mean(f1distpdf(y,nu).*f1distpdf(x(ii)^2-y,nu))*x(ii)^2;
    qq(ii)=2*x(ii)*pp(ii);
end
sum(qq)