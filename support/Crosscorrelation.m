function [C, DOF] = Crosscorrelation(u,v,maxlag)

% v = v-mean(v);
sigma_u = std(u);
sigma_v = std(v);
C = zeros(size(v));
DOF = length(v) - (0:maxlag)';
for lag=0:maxlag;
    v_shift = circshift(v,-lag,1);
    v2 = v.*v_shift;
%     AC(lag+1) = mean(v2(1:(length(v)-lag))); % Crude estimate
    C(lag+1) = sum(v2(1:(length(v)-lag)))/DOF(lag+1); % Proper accounting
end
C = C(1:(maxlag+1))/sigma2;