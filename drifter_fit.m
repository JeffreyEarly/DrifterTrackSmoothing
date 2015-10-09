%------------------------------------------------------------------------------
%
%     Constrained cbs curve fitting in tension
%     Nick Teanby 30/01/07
%
%------------------------------------------------------------------------------
%
%     Function to smooth a set of unevenly spaced x,y data
%	and output the cubic B spline parameters and covariance.
%
%     Allows constraints to be imposed on the smooth curve
%     by using the method of Lagrange multipliers.
%
%	Constraints [optional] can be on y, dy/dx, or d2y/dx2.
%
%	NB. curvature at ends of curve will be constrained to zero by default.
%
%     Tension is applied using a quadratic spring approximation as explained
%	in Teanby 2007.
%
%	  input
%       -----
%	x	float(n)		x data
%	y	float(n)		y data
%	dy	float(n)		y data errors
%	M	int			number of splines to use
%	gamma	float			tension
%
%	  input [optional]
%	  ----------------
%	x0	float(n0)		x constraints
%	y0	float(n0)		constraints [ y , dy/dx, d2y/dx2 ]
%	ctype	float(n0)		constraint type [0=y, 1=grad, 2=curvature]
%
%       output
%       ------
%	m	float(M)		spline parameters
%	cm	float(M,M)		covariance matrix of spline parameters
%
%------------------------------------------------------------------------------
function [m_x,m_y,Cm_x,Cm_y,X,V,A,J] = drifter_fit(t,x,y,dx,dy,T_decorrelation,M,S,v0,lat0, weight_function)

if (length(t) ~= length(x) || length(t) ~= length(y) )
   disp('The time series are not consistent lengths');
   return;
end

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

% M is given as the number of splines on the interval [t(1), t(end)]. We
% need to add extra splines before an after this, depending on order of the
% spline.
M = M + 2*floor(S/2);

% Compute the Coriolis parameter
Omega = 2*pi/86164;
f0 = 2*Omega*sin(lat0*pi/180);

% number of data points
N = length(y);

% Length of series
% T = t(N)-t(1);

% knot spacing
t_knot = (t(N)-t(1))/(M-S);

% spacing and number of points in the quadrature grid
DT = t_knot/500.;
Q = ceil( (t(N)-t(1))/DT + 1 );

gamma = 1/(v0*v0*Q);
% gamma = DT/(dx(1)*dx(1)*(60*30))

% Rows are the N observations
% Columns are the M splines
X = zeros(N,M);
V = zeros(N,M);
A = zeros(N,M);
J = zeros(N,M);
for i=1:N
    for j=1:M
        t_norm=(t(i)-t(1))/t_knot - (j - 1 - floor(S/2));
        X(i,j)=spline(t_norm);
        V(i,j)=spline_t(t_norm);
        A(i,j)=spline_tt(t_norm);
        J(i,j)=spline_ttt(t_norm);
    end
end
V = V/t_knot;
A = A/(t_knot^2);
J = J/(t_knot^3);

% set up F matrix and h vector for constraints
NC = 2;
Fx_x=zeros(NC,M); % force in x, from x constraint
Fx_y=zeros(NC,M); % force in x, from y constraint
Fy_x=zeros(NC,M); % force in y, from x constraint
Fy_y=zeros(NC,M); % force in y, from y constraint
F=zeros(NC,M);

% constrain ends of curve to have zero curvature
for i=1:M
    t_norm=(t(1)-t(1))/t_knot - (i - 1 - floor(S/2));
    Fx_x(1,i) = spline_ttt(t_norm)/(t_knot^3); % x_ttt
    Fy_x(1,i) = f0*spline_tt(t_norm)/(t_knot^2); % + f0 y_tt
    Fx_y(1,i) = -f0*spline_tt(t_norm)/(t_knot^2); % - f0 x_tt
    Fy_y(1,i) = spline_ttt(t_norm)/(t_knot^3); % + y_ttt
    F(1,i) = spline_ttt(t_norm)/(t_knot^3);
    
    t_norm=(t(N)-t(1))/t_knot - (i - 1 - floor(S/2));
    Fx_x(2,i) = spline_ttt(t_norm)/(t_knot^3); % x_ttt
    Fy_x(2,i) = f0*spline_tt(t_norm)/(t_knot^2); % + f0 y_tt
    Fx_y(2,i) = -f0*spline_tt(t_norm)/(t_knot^2); % - f0 x_tt
    Fy_y(2,i) = spline_ttt(t_norm)/(t_knot^3); % + y_ttt
    F(2,i) = spline_ttt(t_norm)/(t_knot^3);
end
h(1) = 0.;
h(2) = 0.;

F_force = 1e9*[Fx_x', Fx_y';Fy_x', Fy_y'];
F_constraint = 1e9*[Fx_x,Fy_x;Fx_y,Fy_y];

% unit test, these should be equivalent to above when lat=0
F_force_known = [F',zeros(M,NC); zeros(M,NC),F'];
F_constraint_known = [F,zeros(NC,M); zeros(NC,M),F];

spline2 = @(t,t1,t2) spline_t(t-t1).*spline_t(t-t2);
spline_spline_t = @(t,t1,t2) spline(t-t1).*spline_t(t-t2);
VV = zeros(M,M);
XV = zeros(M,M);
Mmax = M-2*floor(S/2)-1;
k=0;
for i=1:M
    for j=1:M
        if abs(i-j) <= S
            if (i-S>2) && (i+S<Mmax) && (j-S>2) && (j+S<Mmax)
                VV(i,j) = VV(i-1,j-1);
                XV(i,j) = VV(i-1,j-1);
                k=k+1;
            else
                t1 = (i - 1 - floor(S/2));
                t2 = (j - 1 - floor(S/2));
                VV(i,j) = integral(@(t)spline2(t,t1,t2),0,M-2*floor(S/2)-1);
                XV(i,j) = integral(@(t)spline_spline_t(t,t1,t2),0,M-2*floor(S/2)-1);
            end
        end
    end
end
VV = VV/(t_knot*DT);
XV = XV/DT;

W = WeightingFunction(t,T_decorrelation);
Wx = W*diag(1./(dx.^2));
Wy = W*diag(1./(dy.^2));

dbstop if warning
[m_x,m_y,Cm_x,Cm_y] = ComputeSolution( X, XV, VV, F_force, F_constraint, Wx, Wy, gamma, f0, M, NC, x, y, h );

% if nargin == 10
    error_y_previous = dy;
    rel_error = 1.0;
    repeats = 1;
    while (rel_error > 0.01)
        % These are the deviations of the model from the data
        dx1 = X*m_x - x;
        dy1 = X*m_y - y;
        
%         dx1(find(abs(dx1)<1e-12)) = 1e-12;
%         dy1(find(abs(dy1)<1e-12)) = 1e-12;
      dx1 = max(0.00001*ones(size(dx1)),abs(dx1));
      dy1 = max(0.00001*ones(size(dy1)),abs(dy1));
      ds1 = max(abs(dx1),abs(dy1));
        
        dx2 = dx1./weight_function(dx1./dx);
        dy2 = dy1./weight_function(dy1./dx);
        ds2 = ds1./weight_function(ds1./dx);

        Wx = W*diag(1./(dx2.^2));
        Wy = W*diag(1./(dy2.^2));
        Ws = W*diag(1./(ds2.^2));
        
%          dbstop if warning
        [m_x,m_y,Cm_x,Cm_y] = ComputeSolution( X, XV, VV, F_force, F_constraint, Ws, Ws, gamma, f0, M, NC, x, y, h );

        rel_error = max((dx2-error_y_previous)./dx2);
        error_y_previous=dx2;
        repeats = repeats+1;

        if (repeats == 100)
           disp('Failed to converge after 100 iterations.');
           break;
        end
    end
% end




end

function [m_x,m_y,Cm_x,Cm_y] = ComputeSolution( X, XV, D2, F_force, F_constraint, Wx, Wy, gamma, f0, M, NC, x, y, h )
    % A (and D) matrix:
    % Rows are the N observations
    % Columns are the M splines
    % F is NCxM
    % H is NCx1
    
    % set up inverse theory matrices
%     E_x = A'*Wx*A + gamma*(Dq'*Dq); % MxM
%     E_y = A'*Wy*A + gamma*(Dq'*Dq); % MxM
    E_x = X'*Wx*X + gamma*D2; % MxM
    E_y = X'*Wy*X + gamma*D2; % MxM
    C_x = -f0*gamma*(XV'); % MxM
    C_y = -f0*gamma*(XV); % MxM

    G1 = [[E_x,C_x;C_y,E_y], F_force; F_constraint, zeros(2*NC,2*NC)];
    
%     G1 = [E_x,C_x,Fx';
%           C_y,E_y,Fy';
%           Fx,zeros(NC,NC),zeros(NC,NC);
%           Fy,zeros(NC,NC),zeros(NC,NC)];
    G2 = [X'*Wx*x;X'*Wy*y;h';h'];
    m = G1\G2;
    m_x = m(1:M);
    m_y = m(M+1:2*M);
    
    % model parameter covariance matrix
%     Cm_x = inv(E_x) - inv(E_x)*F'*inv(F*inv(E_x)*F')*F*inv(E_x);
%     Cm_y = inv(E_y) - inv(E_y)*F'*inv(F*inv(E_y)*F')*F*inv(E_y);
Cm_x = zeros(size(E_x));
Cm_y = zeros(size(E_x));
end

function [W] = WeightingFunction(t, T_decorrelation)
    type = 'exponential';
    
    N = length(t);
    
    if (T_decorrelation == 0)
        W = eye(N);
        return;
    end
    
    W=zeros(N,N);
    T=zeros(N,N);
    for i=1:N
        for j=1:N
            W(i,j)=abs(t(i)-t(j));
            T(i,j)= min(min( abs(t(min(i,j))-t(1)), abs(t(max(i,j))-t(end))), T_decorrelation);
        end
    end
    if strcmp(type,'gaussian') == 1
        %T_g2 = -T_decorrelation*T_decorrelation/log(0.01);
        T_g2 = -T.*T/log(0.01);
        W=exp(-W.*W./T_g2);
        W(find(isnan(W)))=1;
    elseif strcmp(type,'exponential') == 1
        T_e = -T_decorrelation/log(0.01);
        W=exp(-W/T_e);
    else
        W=eye(N);
    end
end

% G1_x = [E_x,F';F,zeros(NC:NC)];
% G1_y = [E_y,F';F,zeros(NC:NC)];
% G2_x = [A'*Wx*x;h'];
% G2_y = [A'*Wy*y;h'];
% 
% % invert
% mm_x=inv(G1_x)*G2_x;
% mm_y=inv(G1_y)*G2_y;
% % spline coeffs are first M params
% m_x=mm_x(1:M);
% m_y=mm_y(1:M);


% set up D matrix
% Xq = zeros(Q,M); % Same as the A matrix above, but on the quadrature (q) grid.
% Vq = zeros(Q,M);
% for q=1:Q
%     tq = t(1) + (q-1)*DT;
%     for j=1:M
%         t_norm=(tq-t(1))/t_knot - (j - 1 - floor(S/2));
%         Xq(q,j)=spline(t_norm);
%         Vq(q,j) = spline_t(t_norm)/t_knot;
%     end
% end

