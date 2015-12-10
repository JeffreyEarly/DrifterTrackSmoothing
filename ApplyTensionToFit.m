clear
drifters = load('sample_data/projected_ungridded_rho1_drifters.mat');

 distribution = 'student-t';
%distribution = 'gaussian';

iDrifter = 1;
sigma_gps = 9; % error in meters

S = 3; % order of the spline

if strcmp(distribution,'gaussian')
    p = @(z,sigma) exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi));
    w = @(z)(sigma_gps*sigma_gps);
    gamma2 = 1e9*2.^(0:16)'; % range of tensions we want to explore
elseif strcmp(distribution,'student-t')
    sigma_gps = 9;
    nu = 2.015;
    gamma2 = linspace(2e8,8e11,15)';
    
    nu = 2.033;
    sigma_gps = 14.5;
    gamma2 = linspace(2e9,8e10,16)';
    
    p = @(z,sigma) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2));
    w = @(z)((nu/(nu+1))*sigma_gps^2*(1+z.^2/(nu*sigma_gps^2)));
%     gamma2 = 1e8*2.^(1:16)'; % range of tensions we want to explore
%     gamma2 = linspace(2e8,8e11,15)';
% gamma2 = [121333333333.333];
% gamma2 = [24800000000.0000];
else
    disp('Not a valid distribution')
    return;
end

lat0 = drifters.lat0;
x = drifters.x{iDrifter};
y = drifters.y{iDrifter};
t = drifters.t{iDrifter};

dx = ones(size(x))*sigma_gps;
dy = ones(size(y))*sigma_gps;

t_knot = t;

maxlag = 30;
n = length(x);
coeff=n*(n+2)./(n-(1:maxlag));
variance = zeros(size(gamma2));
ACx = zeros(length(gamma2),maxlag+1);
ACy = zeros(length(gamma2),maxlag+1);
epsilon = zeros(length(gamma2),2*n);
for i=1:length(gamma2)
    [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = drifter_fit_bspline(t,x,y,dx,dy,S,t_knot,[0;gamma2(i)],w);
    X = squeeze(B(:,:,1));
    error_x = X*m_x - x;
    error_y = X*m_y - y;
    epsilon(i,:) = [error_x;error_y];
    variance(i) = mean((epsilon(i,:)).^2);
    ACx(i,:) = Autocorrelation(error_x,maxlag);
    ACy(i,:) = Autocorrelation(error_y,maxlag);
end

% https://en.wikipedia.org/wiki/Ljung?Box_test
% rows test, columns are degrees of freedom
Q = cumsum(ACx(:,2:end).*ACx(:,2:end).*repmat(coeff,length(gamma2),1),2);

Xq = squeeze(Bq(:,:,1));
x_fit = Xq*m_x;
y_fit = Xq*m_y;

figure
subplot(2,2,[1 3])
s = 1/1000;
plot(s*x,s*y), hold on
plot(s*x_fit,s*y_fit,'g')
scatter(s*x,s*y,5)
xlabel('x (km)')
ylabel('y (km)')

subplot(2,2,2)
plot(drifters.t{iDrifter}/3600,s*x), hold on
plot(tq/3600,s*x_fit,'g')
scatter(drifters.t{iDrifter}/3600,s*x,5)
xlabel('t (hours)')
ylabel('x (km)')

subplot(2,2,4)
plot(drifters.t{iDrifter}/3600,s*y), hold on
plot(tq/3600,s*y_fit,'g')
scatter(drifters.t{iDrifter}/3600,s*y,5)
xlabel('t (hours)')
ylabel('y (km)')


var_x = mean((epsilon).^2,2);
sigma_out = sqrt(var_x*(nu-2)/nu );
nu_x_out = 2*var_x./(var_x - sigma_gps*sigma_gps);

figure
xi = (-0500:0.01:0500)';
for i=1:length(gamma2)
    subplot(4,4,i)
    denhist(epsilon(i,:), 500,'b'); hold on
    plot(xi,p(xi,sigma_gps),'LineWidth',2)
    xlim([-50 50])
end

P = [1	1e-7 1e-5	0.001	0.004	0.016	2.706	3.841	5.024	6.635	7.879
2	0.010	0.020	0.051	0.103	0.211	4.605	5.991	7.378	9.210	10.597
3	0.072	0.115	0.216	0.352	0.584	6.251	7.815	9.348	11.345	12.838
4	0.207	0.297	0.484	0.711	1.064	7.779	9.488	11.143	13.277	14.860
5	0.412	0.554	0.831	1.145	1.610	9.236	11.070	12.833	15.086	16.750
6	0.676	0.872	1.237	1.635	2.204	10.645	12.592	14.449	16.812	18.548
7	0.989	1.239	1.690	2.167	2.833	12.017	14.067	16.013	18.475	20.278
8	1.344	1.646	2.180	2.733	3.490	13.362	15.507	17.535	20.090	21.955
9	1.735	2.088	2.700	3.325	4.168	14.684	16.919	19.023	21.666	23.589
10	2.156	2.558	3.247	3.940	4.865	15.987	18.307	20.483	23.209	25.188
11	2.603	3.053	3.816	4.575	5.578	17.275	19.675	21.920	24.725	26.757
12	3.074	3.571	4.404	5.226	6.304	18.549	21.026	23.337	26.217	28.300
13	3.565	4.107	5.009	5.892	7.042	19.812	22.362	24.736	27.688	29.819
14	4.075	4.660	5.629	6.571	7.790	21.064	23.685	26.119	29.141	31.319
15	4.601	5.229	6.262	7.261	8.547	22.307	24.996	27.488	30.578	32.801
16	5.142	5.812	6.908	7.962	9.312	23.542	26.296	28.845	32.000	34.267
17	5.697	6.408	7.564	8.672	10.085	24.769	27.587	30.191	33.409	35.718
18	6.265	7.015	8.231	9.390	10.865	25.989	28.869	31.526	34.805	37.156
19	6.844	7.633	8.907	10.117	11.651	27.204	30.144	32.852	36.191	38.582
20	7.434	8.260	9.591	10.851	12.443	28.412	31.410	34.170	37.566	39.997
21	8.034	8.897	10.283	11.591	13.240	29.615	32.671	35.479	38.932	41.401
22	8.643	9.542	10.982	12.338	14.041	30.813	33.924	36.781	40.289	42.796
23	9.260	10.196	11.689	13.091	14.848	32.007	35.172	38.076	41.638	44.181
24	9.886	10.856	12.401	13.848	15.659	33.196	36.415	39.364	42.980	45.559
25	10.520	11.524	13.120	14.611	16.473	34.382	37.652	40.646	44.314	46.928
26	11.160	12.198	13.844	15.379	17.292	35.563	38.885	41.923	45.642	48.290
27	11.808	12.879	14.573	16.151	18.114	36.741	40.113	43.195	46.963	49.645
28	12.461	13.565	15.308	16.928	18.939	37.916	41.337	44.461	48.278	50.993
29	13.121	14.256	16.047	17.708	19.768	39.087	42.557	45.722	49.588	52.336
30	13.787	14.953	16.791	18.493	20.599	40.256	43.773	46.979	50.892	53.672];
p = [1.0 0.995	0.99	0.975	0.95	0.90	0.10	0.05	0.025	0.01	0.005];
df = P(:,1);
chi2 = [zeros(30,1),P(:,2:end)];

% Compute the probability at each degree of freedom
pQ = zeros(size(Q));
for i=1:size(Q,2)
    pQ(:,i) = interp1(chi2(i,:),p,(Q(:,i)),'linear',0);
end

