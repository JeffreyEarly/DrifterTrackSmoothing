clear
drifters = load('sample_data/projected_ungridded_rho1_drifters.mat');
numDrifters = length(drifters.x);

SplineFactor = 1.5; % Number of data points for each spline
sigma_gps = 10; % error in meters
S = 5; % order of the spline
u_rms = 0.01; % assumed rms velocity of the solution
T_decorrelation = 3*60*60; % forcing decorrelation time
lat0 = drifters.lat0;
lat0 = 0;

t = (0:30*60:drifters.maxExperimentLength)';
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

x = zeros(length(t),numDrifters);
y = zeros(length(t),numDrifters);
u = zeros(length(t),numDrifters);
v = zeros(length(t),numDrifters);

for iDrifter=1:length(drifters.x)
    N = length(drifters.x{iDrifter});
    M = floor(N/SplineFactor); % Number of splines
    dx = ones(size(drifters.x{iDrifter}))*sigma_gps;
    dy = ones(size(drifters.x{iDrifter}))*sigma_gps;
    [mx,my,Cmx,Cmy,X1,V1] = drifter_fit(drifters.t{iDrifter},drifters.x{iDrifter},drifters.y{iDrifter},dx,dy,T_decorrelation,M,S,u_rms,lat0, @(z)(z./(1+0.5*z.*z)));
    
    % Now we create the basis at the desired collocation points.
    M_norm = M+2*floor(S/2);
    t_knot = (drifters.t{iDrifter}(end)-drifters.t{iDrifter}(1))/(M_norm-S);
    X = zeros(Nt,M_norm);
    V = zeros(Nt,M_norm);
    for i=1:Nt
        for j=1:M_norm
            t_norm=(t(i)-drifters.t{iDrifter}(1))/t_knot - (j - 1 - floor(S/2));
            X(i,j)=spline(t_norm);
            V(i,j)=spline_t(t_norm);
        end
    end
    V = V/t_knot;
    
    x(:,iDrifter) = X*mx;
    y(:,iDrifter) = X*my;
    u(:,iDrifter) = V*mx;
    v(:,iDrifter) = V*my;
end

goodIndices = 3:length(t);
t=t(goodIndices,:);
x=x(goodIndices,:);
y=y(goodIndices,:);
u=u(goodIndices,:);
v=v(goodIndices,:);
t = t - t(1);
lon0 = drifters.lat0;
f0 = drifters.f0;

save('griddedRho2DriftersWithTensionFits.mat', 't', 'x', 'y', 'lat0', 'lon0', 'f0')
