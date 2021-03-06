clear
drifters = load('sample_data/projected_ungridded_rho1_drifters.mat');
lat0 = drifters.lat0;
Omega = 2*pi/86164;
f0 = 2*Omega*sin(lat0*pi/180);

 distribution = 'student-t';
 %distribution = 'gaussian';

iDrifter = 8;
sigma = 9; % error in meters

S = 4; % order of the spline
K = S+1;

if strcmp(distribution,'gaussian')
    p = @(z) exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi));
    w = @(z)(sigma*sigma);
    gamma2 = 1e9*2.^(0:15)'; % range of tensions we want to explore
    gamma2 = 10.^(linspace(log10(2.5e2),log10(2.5e8),16)');
elseif strcmp(distribution,'student-t')
     nu = 2.005; sigma = 8;
%    nu = 2.019; sigma = 15;
    var = 158.^2;
    sigma = 20;
    nu = 2*var/(var-sigma*sigma)
    a_rms =  10.^(linspace(log10(2e-6),log10(1e-4),16)'); % optimal is 1.4e-5 m/s^2
    
      a_rms =  10.^(linspace(log10(5e-10),log10(1e-7),16)'); % optimal is 7e-9 m/s^3
     
%     a_rms =  10.^(linspace(log10(5e-14),log10(1e-10),16)'); % optimal is 3.8e-12 m/s^4
%     
    p = @(z) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2));
    w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));
    rho = @(z) -log(sqrt(2)*gamma((nu+1)/2)./(sqrt(nu)*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2)));
else
    disp('Not a valid distribution')
    return;
end

Nall = 0;
for iDrifter = 1:9
    Nall = Nall + length(drifters.x{iDrifter});
end


maxlag = 30;


ACx = zeros(length(a_rms),maxlag+1);
ACy = zeros(length(a_rms),maxlag+1);
epsilon = zeros(length(a_rms),2*Nall);
Ndrifters = length(drifters.x);
twiddle_factors = zeros(size(length(a_rms),Ndrifters));
a_rms_out = zeros(size(length(a_rms),Ndrifters));
for i=1:length(a_rms)
    fprintf('tension %d (of %d)\n',i,length(a_rms));
    iEpsilon = 1;
    for iDrifter = 1:Ndrifters
        x = drifters.x{iDrifter};
        y = drifters.y{iDrifter};
        t = drifters.t{iDrifter};
        N = length(t);
        t_knot = t;
        gamma_in = N./(a_rms(i).*a_rms(i)*(t(end)-t(1)));
        [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = drifter_fit_bspline(t,x,y,ones(size(x))*sigma,ones(size(y))*sigma,S,t_knot,gamma_in,w);
        X = squeeze(B(:,:,1));
        error_x = X*m_x - x;
        error_y = X*m_y - y;
        epsilon(i,iEpsilon:(iEpsilon+2*N-1)) = [error_x;error_y];
        iEpsilon = iEpsilon + 2*N + 1;
        ACx(i,:) = ACx(i,:) + Autocorrelation(error_x,maxlag)';
        ACy(i,:) = ACy(i,:) + Autocorrelation(error_y,maxlag)';
        
        jx = squeeze(Bq(:,:,S))*m_x;
        jy = squeeze(Bq(:,:,S))*m_y;
        j = jx.*jx+jy.*jy;
        twiddle_factors(i,iDrifter) = gamma_in/(2*N/(sum(j)*(tq(2)-tq(1))));
        a_rms_out(i,iDrifter) = sqrt( sum(j)*(tq(2)-tq(1))/(tq(end)-tq(1)) );
    end
end

[a_rms, sqrt(mean(a_rms_out.*a_rms_out,2))]

if strcmp(distribution,'student-t')
    var_x = mean((epsilon).^2,2);
    sigma_out = sqrt(var_x*(nu-2)/nu );
    nu_x_out = 2*var_x./(var_x - sigma*sigma);
end


ACx = ACx/Ndrifters;
ACy = ACy/Ndrifters;
AC = (ACx+ACy)/2;

histwidth = 100;
nbins = 50;
binwidth = histwidth/nbins;
edges = [-3000;((-histwidth/2+binwidth):binwidth:(histwidth/2-binwidth))';3000]; % edges of all bins
binleft = linspace((-histwidth/2),(histwidth/2-binwidth),nbins)'; % left edges of all the ones we'll show

edgedensity = integral(p,(histwidth/2-binwidth),3000)/binwidth;

figure
xi_left = linspace(-histwidth/2,-histwidth/2+binwidth,10)';
xi_mid = linspace(-histwidth/2+binwidth,histwidth/2-binwidth,100)';
xi_mid = xi_mid(2:(end-1));
xi_right = linspace(histwidth/2-binwidth,histwidth/2,10)';
xi = [xi_left;xi_mid;xi_right];
denfunc = edgedensity*ones(size(xi));
denfunc(11:108) = p(xi_mid);
for i=1:length(a_rms)
    subplot(4,4,i)
    %     denhist(epsilon(i,:), 500,'b'); hold on
    
    count = histcounts(epsilon(i,:),edges)';
    bar(binleft, count/(length(epsilon(i,:))*binwidth), 'histc'); hold on;
    
    plot(xi,denfunc,'LineWidth',2,'Color','magenta')
    xlim([-50 50])
end

% Confirmation that the we're summing to 1.0
% trapz(xi,denfunc);
% trapz(binleft, count/(length(epsilon(i,:))*binwidth));

% https://en.wikipedia.org/wiki/Ljung?Box_test
% rows test, columns are degrees of freedom
n = length(x);
coeff=n*(n+2)./(n-(1:maxlag));
AC = (ACx+ACy)/2;
Qx = cumsum(ACx(:,2:end).*ACx(:,2:end).*repmat(coeff,length(a_rms),1),2);
Qy = cumsum(ACy(:,2:end).*ACy(:,2:end).*repmat(coeff,length(a_rms),1),2);
Q = cumsum(AC(:,2:end).*AC(:,2:end).*repmat(coeff,length(a_rms),1),2);

chicheck = sum((epsilon/sigma).^2,2)/length(t);

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
prob = [1.0 0.995	0.99	0.975	0.95	0.90	0.10	0.05	0.025	0.01	0.005];
df = P(:,1);
chi2 = [zeros(30,1),P(:,2:end)];

figure
plot(Q')
hold on, plot(P(:,6)','LineWidth',2)

% subplot(1,2,1)
% plot(Qx')
% hold on, plot(P(:,6)','LineWidth',2)
% 
% subplot(1,2,2)
% plot(Qy')
% hold on, plot(P(:,6)','LineWidth',2)


% Compute the probability at each degree of freedom
pQ = zeros(size(Q));
for i=1:size(Q,2)
    pQ(:,i) = interp1(chi2(i,:),prob,(Q(:,i)),'linear',0);
end
