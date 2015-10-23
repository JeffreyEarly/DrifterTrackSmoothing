clear
drifters = load('sample_data/projected_ungridded_rho1_drifters.mat');

iDrifter = 1;
SplineFactor = 1.1; % Number of data points for each spline
sigma_gps = 10; % error in meters
S = 5; % order of the spline
u_rms = 0.15; % assumed rms velocity of the solution
T_decorrelation = 0; %2.5*60*60; % forcing decorrelation time

% Note that big decorrelation times are causing the the solvability issue.
% for some reason the the x3 points are not collocating.

lat0 = drifters.lat0;
x = drifters.x{iDrifter};
y = drifters.y{iDrifter};
t = drifters.t{iDrifter};

error_y = ones(size(y))*sigma_gps; % This is the error, in kilometers, of the GPS fits.
N = length(y);
M = floor(N/SplineFactor); % Number of splines

dx = ones(size(x))*sigma_gps;
dy = ones(size(y))*sigma_gps;

[mx,my,Cmx,Cmy,X,V] = forcing_fit_cauchy(t,x,y,dx,dy,T_decorrelation,M,S,1000,0,@(z)(z./(1+0.5*z.*z)));
x1 = X*mx;
y1 = X*my;
u1 = V*mx;
v1 = V*my;

% Now we create the basis at the desired collocation points.
t = (0:5*60:drifters.maxExperimentLength)';
Nt = length(t);
if S == 3
    addpath('./cubic_splines');
    spline = @(t) cspline(t);
    spline_t = @(t) cspline_t(t);
    spline_tt = @(t) cspline_tt(t);
    spline_ttt = @(t) cspline_ttt(t);
elseif S == 5
    addpath('./quintic_splines');
    spline = @(t) qspline(t);
    spline_t = @(t) qspline_t(t);
    spline_tt = @(t) qspline_tt(t);
    spline_ttt = @(t) qspline_ttt(t);
else
    disp('Whoops! I only know how to deal with splines of order 3 and 5.')
    return;
end
M_norm = M+2*floor(S/2);
t_knot = (drifters.t{iDrifter}(end)-drifters.t{iDrifter}(1))/(M_norm-S);
X = zeros(Nt,M_norm);
V = zeros(Nt,M_norm);
A = zeros(Nt,M_norm);
J = zeros(Nt,M_norm);
for i=1:Nt
    for j=1:M_norm
        t_norm=(t(i)-drifters.t{iDrifter}(1))/t_knot - (j - 1 - floor(S/2));
        X(i,j)=spline(t_norm);
        V(i,j)=spline_t(t_norm);
        A(i,j)=spline_tt(t_norm);
        J(i,j)=spline_ttt(t_norm);
    end
end
V = V/t_knot;
A = A/t_knot^2;
J = J/t_knot^3;

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

return
Diff2 = FiniteDifferenceMatrix(2,drifters.t{iDrifter},2,2,10);
Diff2 = Diff2(2:(N-1),:);
t_obs_force = drifters.t{iDrifter}(2:(end-1));
figure
subplot(2,1,1)
scatter(t_obs_force,Diff2*x), hold on
plot(t,A*mx)
subplot(2,1,2)
scatter(t_obs_force,Diff2*y), hold on
plot(t,A*my)

