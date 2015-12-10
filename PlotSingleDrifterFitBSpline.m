clear
drifters = load('sample_data/projected_ungridded_rho1_drifters.mat');

iDrifter = 1;
SplineFactor = 0.5; % Number of data points for each spline
sigma_gps = 9; % error in meters
nu = 2.015;
S = 5; % order of the spline
u_rms = 1e-13; % assumed rms velocity of the solution
T_decorrelation = 0; %2.5*60*60; % forcing decorrelation time

% SplineFactor = 1.5; % Number of data points for each spline
% sigma_gps = 9; % error in meters
% nu = 2.015;
% u_rms = 1; % assumed rms velocity of the solution

iDrifter = 1;
SplineFactor = 5; % Number of data points for each spline
sigma_gps = 9; % error in meters
nu = 2.015;
S = 4; % order of the spline
u_rms = 1e-13; % assumed rms velocity of the solution
T_decorrelation = 0; %2.5*60*60; % forcing decorrelation time


% sigma=8.5 and nu = 2.01

lat0 = drifters.lat0;
x = drifters.x{iDrifter};
y = drifters.y{iDrifter};
t = drifters.t{iDrifter};

N = length(y);
M = floor(N/SplineFactor); % Number of splines

dx = ones(size(x))*sigma_gps;
dy = ones(size(y))*sigma_gps;

p_g = @(z,sigma) exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi));
p_t = @(z,sigma) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2));

% Gaussian and t-distribution weighting functions.
w_g = @(z)(sigma_gps*sigma_gps);
w_t = @(z)((nu/(nu+1))*sigma_gps^2*(1+z.^2/(nu*sigma_gps^2)));


t_knot = t(1:1:end);
[m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = drifter_fit_bspline(t,x,y,dx,dy,S,t_knot,[0,0,0],w_g);
X = squeeze(B(:,:,1));
V = squeeze(B(:,:,2));
A = squeeze(B(:,:,3));
x1 = X*m_x;
y1 = X*m_y;
u1 = V*m_x;
v1 = V*m_y;
ax = A*m_x;
ay = A*m_y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Error distribution plots
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure

 % position error
epsilon = [x-x1;y-y1];
var_x = mean((epsilon).^2);
sigma_out = sqrt(var_x*(nu-2)/nu );
nu_x_out = 2*var_x/(var_x - sigma_gps*sigma_gps);
xi = (-0500:0.01:0500)';

subplot(1,3,1)
denhist(epsilon, 500,'b'); hold on
plot(xi,p_g(xi,sigma_gps),'LineWidth',2)
plot(xi,p_t(xi,sigma_gps),'LineWidth',2)
vlines([-sigma_out sigma_out],'r--')
xlim([-250 250])
xlabel('model error (meters)')
title(sprintf('\\sigma_{in}=%1.2f, \\nu_{in}=%1.3f\n \\sigma_{out}=%1.2f, \\nu_{out}=%1.4g',sigma_gps,nu, sigma_out,nu_x_out))


% velocity error
[Diff1,t_v] = FiniteDifferenceMatrixNoBoundary(1, t, 1);
epsilon_v = [Diff1*x-Diff1*x1;Diff1*y-Diff1*x1];
xi = (-0.20:0.001:0.20)';
dt = sqrt(mean(diff(t).^2));
sigma_v = sqrt(2)*sigma_gps/dt;
var_v = mean((epsilon_v).^2);
sigma_v_out = (dt/sqrt(2))*sqrt(var_v*(nu-2)/nu);
nu_v_out = 2*var_v/(var_v - sigma_v*sigma_v);

subplot(1,3,2)
denhist(epsilon_v, 500,'b'); hold on
plot(xi,p_g(xi,sigma_v),'LineWidth',2)
plot(xi,p_t(xi,sigma_v),'LineWidth',2)
vlines([-sigma_v_out sigma_v_out],'r--')
xlim([-0.2 0.2])
xlabel('model error (meters/second)')
title(sprintf('\\sigma_{out}=%1.2f, \\nu_{out}=%1.4g', sigma_v_out, nu_v_out))




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

return;

% acceleration error
[Diff2,t_a] = FiniteDifferenceMatrixNoBoundary(2, t, 1);

epsilon_a = [Diff2*x-ax;Diff2*y-ay];
xi = (-5e-4:1e-6:5e-4)';
sigma_a = sqrt(6)*sigma_gps/dt^2;
var_a = mean((epsilon_a).^2);
sigma_a_out = (dt^2/sqrt(6))*sqrt(var_a*(nu-2)/nu);
nu_a_out = 2*var_a/(var_a - sigma_a*sigma_a);

subplot(1,3,3)
denhist(epsilon_a, 500,'b'); hold on
plot(xi,p_g(xi,sigma_a),'LineWidth',2)
plot(xi,p_t(xi,sigma_a),'LineWidth',2)
xlim([-2e-4 2e-4])
xlabel('model error (meters/second^2)')
title(sprintf('\\sigma_{out}=%1.2f, \\nu_{out}=%1.4g', sigma_a_out, nu_a_out))

[Diff3,t_j] = FiniteDifferenceMatrixNoBoundary(3, t, 1);
sigma_j = sqrt(4*3*2)*sigma_gps/dt^3;

% Now we create the basis at the desired collocation points.
t = (0:5*60:drifters.maxExperimentLength)';
Nt = length(t);
t_norm = zeros(Nt,M_norm);
for i=1:Nt
    for j=1:M_norm
        t_norm(i,j)=(t(i)-drifters.t{iDrifter}(1))/t_knot - (j - 1 - floor(S/2));
    end
end
X=spline(t_norm);
V = spline_t(t_norm)/t_knot;
A = spline_tt(t_norm)/t_knot^2;
J = spline_ttt(t_norm)/t_knot^3;


figure
subplot(2,2,[1 3])
s = 1/1000;
plot(s*x,s*y), hold on
plot(s*X*mx,s*X*my,'g')

scatter(s*x,s*y,5)
% xlim([-500 7500])
% ylim([2000 10000])
xlabel('x (km)')
ylabel('y (km)')
title('Fit with Cauchy errors, no tension')

subplot(2,2,2)
plot(drifters.t{iDrifter}/3600,s*x), hold on
plot(t/3600,s*X*mx,'g')
scatter(drifters.t{iDrifter}/3600,s*x,5)
% xlim([8.295e6 1.01e7])
xlabel('t (hours)')
ylabel('x (km)')

subplot(2,2,4)
plot(drifters.t{iDrifter}/3600,s*y), hold on
plot(t/3600,s*X*my,'g')
scatter(drifters.t{iDrifter}/3600,s*y,5)
% xlim([8.295e6 1.01e7])
xlabel('t (hours)')
ylabel('y (km)')

stdColor = [0.8 0.8 0.8];
stdEdgeColor = 'none';
% t_obs=drifters.t{iDrifter};
% figure
% subplot(2,1,1)
% error_plot = [x+sigma_gps*ones(size(t_obs)); flip(x-sigma_gps*ones(size(t_obs)),1)];
% t_error = [t_obs; flip(t_obs,1)];
% fill(t_error/3600,error_plot*s,stdColor,'EdgeColor',stdEdgeColor), hold on
% scatter(drifters.t{iDrifter}/3600,s*x,5, 'k', 'filled')

figure
subplot(2,1,1)
scatter(t_v/3600,Diff1*x), hold on
errorbar(t_v/3600,Diff1*x,sigma_v*ones(size(t_v)))
plot(t/3600,V*mx)
subplot(2,1,2)
scatter(t_v/3600,Diff1*y), hold on
errorbar(t_v/3600,Diff1*y,sigma_v*ones(size(t_v)))
plot(t/3600,V*my)


Sigma_x = sigma_v*ones(size(t_v));
Sigma_y = sigma_v*ones(size(t_v));
for i=1:length(Sigma_x)
    Sigma_x(i) = sum(Diff1(i,:).*Diff1(i,:).*dx2'.*dx2');
    Sigma_y(i) = sum(Diff1(i,:).*Diff1(i,:).*dy2'.*dy2');
end
Sigma_x = sqrt(Sigma_x);
Sigma_y = sqrt(Sigma_y);

figure
subplot(2,1,1)
error_plot = [Diff1*x+Sigma_x; flip(Diff1*x-Sigma_x,1)];
t_v_error = [t_v; flip(t_v,1)];
fill(t_v_error/3600,error_plot,stdColor,'EdgeColor',stdEdgeColor), hold on
plot(t/3600,V*mx,'cyan','LineWidth',2)
scatter(t_v/3600,Diff1*x,5, 'k', 'filled')
ylim([-0.3 0.3])


figure
subplot(2,1,1)
scatter(t_a/3600,Diff2*x), hold on
errorbar(t_a/3600,Diff2*x,sigma_a*ones(size(t_a)))
plot(t/3600,A*mx)
subplot(2,1,2)
scatter(t_a/3600,Diff2*y), hold on
errorbar(t_a/3600,Diff2*y,sigma_a*ones(size(t_a)))
plot(t/3600,A*my)


figure
subplot(2,2,1)
denhist(diff(x),500,'b');
hold on
vlines([-sigma_gps sigma_gps],'r--');
xlim([-10*sigma_gps 10*sigma_gps])
title(sprintf('%2.2f%%',sum(abs(diff(x))>sigma_gps)/length(diff(x))))

% These are points that are uncorrelated
D = zeros(N-3,N-1);
for i=1:size(D,1)
   D(i,i) = -1;
   D(i,i+2) = 1;
end
DD1x = D*Diff1*x;
subplot(2,2,2)
denhist(DD1x,500,'b');
hold on
vlines([-sigma_v sigma_v],'r--');
xlim([-10*sigma_v 10*sigma_v])
title(sprintf('%2.2f%%',sum(abs(DD1x)>sigma_v)/length(DD1x)))

% These are points that are uncorrelated
D = zeros(N-5,N-2);
for i=1:size(D,1)
   D(i,i) = -1;
   D(i,i+3) = 1;
end
DD2x = D*Diff2*x;
subplot(2,2,3)
denhist(DD2x,500,'b');
hold on
vlines([-sigma_a sigma_a],'r--');
xlim([-10*sigma_a 10*sigma_a])
title(sprintf('%2.2f%%',sum(abs(DD2x)>sigma_a)/length(DD2x)))


% These are points that are uncorrelated
D = zeros(N-7,N-3);
for i=1:size(D,1)
   D(i,i) = -1;
   D(i,i+4) = 1;
end
DD3x = D*Diff3*x;
subplot(2,2,4)
denhist(DD3x,500,'b');
hold on
vlines([-sigma_j sigma_j],'r--');
xlim([-10*sigma_j 10*sigma_j])
title(sprintf('%2.2f%%',sum(abs(DD3x)>sigma_j)/length(DD3x)))

% Diff1 = FiniteDifferenceMatrix(1,drifters.t{iDrifter},1,1,1);
% Diff2 = FiniteDifferenceMatrix(2,drifters.t{iDrifter},2,2,2);
% Omega = 2*pi/86164;
% f0 = 2*Omega*sin(lat0*pi/180);
% f_x = A*mx + f0*V*my;
% f_y = A*my - f0*V*mx;
% F_x = Diff2*x + f0*Diff1*y;
% F_y = Diff2*y - f0*Diff1*x;
% 
% figure
% subplot(2,1,1)
% scatter(drifters.t{iDrifter}/3600,F_x), hold on
% plot(t/3600,f_x)
% subplot(2,1,2)
% scatter(drifters.t{iDrifter}/3600,F_y), hold on
% plot(t/3600,f_y)



