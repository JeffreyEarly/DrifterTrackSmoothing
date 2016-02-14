function [qq] = tdistpdf(r, nu)
% numerically computed PDF of r=\sqrt{x^+y^2}, where x and y are
% independent t(nu) distributions
% r pdf vector
% output vector qq is the pdf at locations given by r

nonzeroIndices = find(r>0);
pp=zeros(size(r)); qq=zeros(size(r));
for j = 1:length(nonzeroIndices)
    ii = nonzeroIndices(j);
    gridX = r(ii)^2/400000;
    y=gridX/2:gridX:r(ii)^2; %-gridX/2;
    pp(ii)=mean(f1distpdf(y,nu).*f1distpdf(r(ii)^2-y,nu))*r(ii)^2;
    qq(ii)=2*r(ii)*pp(ii);
end