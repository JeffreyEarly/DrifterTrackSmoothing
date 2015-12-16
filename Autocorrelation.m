function AC = Autocorrelation(v,maxlag)

% v = v-mean(v);
% sigma2 = std(v).^2;
sigma2 = mean(v.*v);
AC = zeros(size(v));
for lag=0:maxlag;
    v_shift = circshift(v,-lag,1);
    v2 = v.*v_shift;
    AC(lag+1) = mean(v2(1:(length(v)-lag))); % Using the mean here makes this a *biased* estimate
end
AC = AC(1:(maxlag+1))/sigma2;