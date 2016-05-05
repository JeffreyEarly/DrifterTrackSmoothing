function [C, DOF] = Crosscorrelation(u,v,maxlag)

n = min(length(u),length(v));
u=u(1:n);
v=v(1:n);

sigma_u = std(u);
sigma_v = std(v);
C = zeros(size(v));
DOF = n - (0:maxlag)';
for lag=0:maxlag;
    v_shift = circshift(v,-lag,1);
    v2 = u.*v_shift;
%     AC(lag+1) = mean(v2(1:(length(v)-lag))); % Crude estimate
    C(lag+1) = sum(v2(1:(n-lag)))/DOF(lag+1); % Proper accounting
end
C = C(1:(maxlag+1))/(sigma_u*sigma_v);