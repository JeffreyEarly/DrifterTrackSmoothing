%addpath('../support/splines/')
lat0 = 0; % Set to 0 for no Coriolis
Omega = 2*pi/86164;
f0 = 2*Omega*sin(lat0*pi/180);

N = 3; % Number of observations
t=linspace(0,86400/1,N);
x=[0;500;1000];
y=[0;500;0];

% t = [0; 86400/3; 86400];
% x=[0;500;1000];
% y=[0;0;0];

dx=10*ones(size(x));
dy=10*ones(size(y));
a0=20;
M=30; % Number of interior knot points (need two extras for end points)
W=eye(N);
S = 5;
[m_x,m_y,Cm_x,Cm_y,X,V,A,J] = drifter_fit_lagrangian(t,x,y,dx,dy,W,M,S,a0,lat0,@(z)(z));

Nt=51;
M = M+2*floor(S/2);
tt = linspace(0,86400,Nt)';
t_knot = (t(end)-t(1))/(M-S);

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

X = zeros(Nt,M);
V = zeros(Nt,M);
A = zeros(Nt,M);
J = zeros(Nt,M);
for i=1:Nt
    for j=1:M
        t_norm=(tt(i)-t(1))/t_knot - (j - 1 - floor(S/2));
        X(i,j)=spline(t_norm);
        V(i,j)=spline_t(t_norm);
        A(i,j)=spline_tt(t_norm);
        J(i,j)=spline_ttt(t_norm);
    end
end

xi = X*m_x;
eta = X*m_y;

% xi = xi(1:M/2);
% eta = eta(1:M/2);

figure
plot(xi,eta)
hold on
scatter(xi,eta)

figure
subplot(4,1,1)
plot(tt/3600,[X*m_x, X*m_y])
ylabel('distance (m)')
subplot(4,1,2)
plot(tt/3600,[V*m_x, V*m_y])
ylabel('speed (m/s)')
subplot(4,1,3)
plot(tt/3600,[A*m_x, A*m_y])
ylabel('acceleration (m/s^2)')
subplot(4,1,4)
plot(tt/3600,[J*m_x, J*m_y])
ylabel('jerk (m/s^3)')
xlabel('time (hours)')