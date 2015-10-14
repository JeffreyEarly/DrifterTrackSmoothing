N = 100000;
x = randn(N,1);
v = diff(diff(x));

ac = zeros(size(v));
for lag=0:30;
    v_shift = circshift(v,-lag,1);
    v2 = v.*v_shift;
    ac(lag+1) = mean(v2(1:(length(v)-lag)));
end