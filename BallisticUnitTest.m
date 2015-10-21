g = -9.8;
v0 = 500;
u0 = 500;
t_max = -2*v0/g;
% t_max = 100;

N =7; % Number of observations
t=linspace(0,t_max,N)';
x=u0*t;
y=g*t.*t/2 + v0*t;

sigma = 100;
x = x + randn(size(x))*sigma*0;
y = y + randn(size(y))*sigma*0;

dx=sigma*ones(size(x));
dy=sigma*ones(size(y));
a0=2*abs(g);
a0 = 1000;
% a0 = 1;
M=N+0*(N-1); % Number of interior knot points (need two extras for end points)
%M = 6;
W=eye(N);
S = 5;
[m_x,m_y,Cm_x,Cm_y,X,V,A,J] = forcing_fit_cauchy(t,x,y,dx,dy,0,M,S,a0,0,@(z)(z));

Nt=3*M;
M = M+2*floor(S/2);
tt = linspace(0,t_max,Nt)';
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
V = V/t_knot;
A = A/(t_knot^2);
J = J/(t_knot^3);

xi = X*m_x;
eta = X*m_y;

% xi = xi(1:M/2);
% eta = eta(1:M/2);

% PEx = zeros(length(tt),N);
% PEy = zeros(length(tt),N);
% for i=1:N
%     PEx(:,i) = (interp1q(tt,xi,t(i))-x(i))*(xi-x(i));
%     PEy(:,i) = (interp1q(tt,eta,t(i))-y(i))*(eta-y(i));
% end
% figure
% plot(tt,PEx)
% figure
% plot(tt,PEy)

figure
plot(xi,eta)
hold on
scatter(xi,eta)
scatter(x,y,5^2,'k')

figure
subplot(4,1,1)
plot(tt,[X*m_x, X*m_y])
xlim([min(t) max(t)])
ylabel('distance (m)')
subplot(4,1,2)
plot(tt,[V*m_x, V*m_y])
xlim([min(t) max(t)])
ylabel('speed (m/s)')
subplot(4,1,3)
plot(tt,[A*m_x, A*m_y])
hold on
scatter(t,zeros(size(t)))
xlim([min(t) max(t)])
ylabel('acceleration (m/s^2)')
subplot(4,1,4)
plot(tt,[J*m_x, J*m_y])
xlim([min(t) max(t)])
ylabel('jerk (m/s^3)')
xlabel('time (hours)')