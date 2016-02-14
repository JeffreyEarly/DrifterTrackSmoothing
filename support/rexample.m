%% EXAMPLE
nu = 5; delta = .01; mx = 10;
% nu: degrees of freedom parameter
% delta: increments in the desired pdf vector
% mx: maximum x-axis value of the desired pdf vector
MC=100000; % number of Monte Carlo Repeats for simulated histograms
BR=100; % number of histogram bars
xxx = trnd(nu,1,MC); yyy = trnd(nu,1,MC);
rrr=sqrt(xxx.^2+yyy.^2); % here is your samples for r
figure; % plotting the data
[n,xout] = hist(min(mx,rrr),BR);
bar(xout,BR/min(mx,max(rrr))*n/sum(n),'BarWidth', 1); %relative frequency is n/sum(n)
qq=rdistpdf(nu,delta,mx); x=delta:delta:mx; % here is the numerically computed pdf for r
hold on; plot(x,qq,'r') % plotting the computed pdf
hold on; plot(x,raylpdf(x),'g') % rayleigh approximation for gaussian
legend('samples','num. distribution','Rayleigh')