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
function [m_x,m_y,Cm_x,Cm_y,X,V,A,J,dx2,dy2] = forcing_fit_cauchy(t,x,y,dx,dy,T_decorrelation,M,S,v0,lat0, weight_function)

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
    spline_tttt = @(t) qspline_tttt(t);
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

% Rows are the N observations
% Columns are the M splines
t_norm = zeros(N,M);
for i=1:N
    for j=1:M
        t_norm(i,j)=(t(i)-t(1))/t_knot - (j - 1 - floor(S/2));   
    end
end
X=spline(t_norm);
A=spline_tt(t_norm)/(t_knot^2);
J=spline_ttt(t_norm)/(t_knot^3);

% Same as the V matrix above, but on the quadrature (q) grid.
t_norm = zeros(Q,M);
for q=1:Q
    tq = t(1) + (q-1)*DT;
    for j=1:M
        t_norm(q,j)=(tq-t(1))/t_knot - (j - 1 - floor(S/2));
    end
end
Vq=spline_t(t_norm)/t_knot;
Jq=spline_ttt(t_norm)/t_knot^3;
Jounce=spline_tttt(t_norm)/t_knot^4;

% Differentiation matrix for velocity, and velocity grid t_v
[Diff1,t_v] = FiniteDifferenceMatrixNoBoundary(1, t, 1);
N_v = length(t_v);
tv_norm = zeros(N_v,M);
for i=1:N_v
    for j=1:M
        tv_norm(i,j)=(t_v(i)-t_v(1))/t_knot - (j - 1 - floor(S/2));   
    end
end
V=spline_t(tv_norm)/t_knot;

% Differentiation matrix for acceleration, and acceleration grid t_a
[Diff2,~] = FiniteDifferenceMatrixNoBoundary(2, t, 1);
A = A(2:(N-1),:);

% set up F matrix and h vector for constraints
NC = 2;
F=zeros(NC,M);

% constrain ends of curve to have zero curvature
for i=1:M
    t_norm=(t(1)-t(1))/t_knot - (i - 1 - floor(S/2));
    F(1,i) = spline_ttt(t_norm)/(t_knot^3);
    t_norm=(t(N)-t(1))/t_knot - (i - 1 - floor(S/2));
    F(2,i) = spline_ttt(t_norm)/(t_knot^3);
end
hx(1) = 0;
hx(2) = 0;
hy(1) = 0;
hy(2) = 0;

% 2nd order 1st derivative matrix, no boundary terms
% Diff1b = FiniteDifferenceMatrix(1,t,1,1,10);
% Diff1b = Diff1b(2:(N-1),:);

% 2nd order 2nd derivative matrix, no boundary terms
% Diff2b = FiniteDifferenceMatrix(2,t,2,2,10);
% Diff2b = Diff2b(2:(N-1),:);

Sigma_x=zeros(N-1,N-1);
Sigma_y=zeros(N-1,N-1);
for i=1:size(Sigma_x,1)
    for j=1:size(Sigma_x,2)
        Sigma_x(i,j) = sum(Diff1(i,:).*Diff1(j,:).*dx'.*dx');
        Sigma_y(i,j) = sum(Diff1(i,:).*Diff1(j,:).*dy'.*dy');
    end
end

Sigma_xx=zeros(N-2,N-2);
Sigma_yy=zeros(N-2,N-2);
for i=1:size(Sigma_xx,1)
    for j=1:size(Sigma_xx,2)
        Sigma_xx(i,j) = sum(Diff2(i,:).*Diff2(j,:).*dx'.*dx');
        Sigma_yy(i,j) = sum(Diff2(i,:).*Diff2(j,:).*dy'.*dy');
    end
end

WXx = diag(1./(dx.^2));
WXy = diag(1./(dy.^2));

WVx = inv(Sigma_x);
WVy = inv(Sigma_y);

WAx = inv(Sigma_xx);
WAy = inv(Sigma_yy);

% Acceleration, but w/o the first and last point
% Vf = V(2:(N-1),:);
% Af = A(2:(N-1),:);

% Tension
T = gamma*(Jq'*Jq);
%T = gamma*(Vq'*Vq);
T = gamma*(Jounce'*Jounce);
x_t = Diff1*x;
y_t = Diff1*y;
x_tt = Diff2*x;
y_tt = Diff2*y;

% [m_x,m_y,Cm_x,Cm_y] = ComputeSolution( X, V, A, T, WXx, WXy, WVx, WVy, WAx, WAy, F, hx, hy, f0, M, NC, x, y, x_t, y_t, x_tt, y_tt );
[m_x,m_y,Cm_x,Cm_y] = ComputeSolutionWithElim( X, V, A, T, WXx, WXy, Sigma_x, Sigma_y, Sigma_xx, Sigma_yy, F, hx, hy, f0, M, NC, x, y, x_t, y_t, x_tt, y_tt );

error_y_previous = dy;
rel_error = 1.0;
repeats = 1;
while (rel_error > 0.01)
    dx1 = X*m_x - x;
    dy1 = X*m_y - y;

    dx2 = weight_function(dx1);
    dy2 = weight_function(dy1);

    for i=1:size(Sigma_x,1)
        for j=1:size(Sigma_x,2)
            Sigma_x(i,j) = sum(Diff1(i,:).*Diff1(j,:).*dx2'.*dx2');
            Sigma_y(i,j) = sum(Diff1(i,:).*Diff1(j,:).*dy2'.*dy2');
        end
    end
    
    for i=1:size(Sigma_xx,1)
        for j=1:size(Sigma_xx,2)
            Sigma_xx(i,j) = sum(Diff2(i,:).*Diff2(j,:).*dx2'.*dx2');
            Sigma_yy(i,j) = sum(Diff2(i,:).*Diff2(j,:).*dy2'.*dy2');
        end
    end

%     WAx = inv(Sigma_x);
%     WAy = inv(Sigma_y);

    WXx = diag(1./(dx2.^2));
    WXy = diag(1./(dy2.^2));
    
    %         dbstop if warning
%     [m_x,m_y,Cm_x,Cm_y] = ComputeSolution( X, V, Af, T, WXx, WXy, WAx, WAy, F, hx, hy, f0, M, NC, x, y, x_t, y_t );
    [m_x,m_y,Cm_x,Cm_y] = ComputeSolutionWithElim( X, V, A, T, WXx, WXy, Sigma_x, Sigma_y, Sigma_xx, Sigma_yy, F, hx, hy, f0, M, NC, x, y, x_t, y_t, x_tt, y_tt );
    
    rel_error = max((dx2-error_y_previous)./dx2);
    error_y_previous=dx2;
    repeats = repeats+1;
    
    if (repeats == 100)
        disp('Failed to converge after 100 iterations.');
        break;
    end
    
end

end


function [m_x,m_y,Cm_x,Cm_y] = ComputeSolution( X, V, A, T, WXx, WXy, WVx, WVy, WAx, WAy, F, hx, hy, f0, M, NC, x, y, x_t, y_t, x_tt, y_tt )

EX_x = X'*WXx*X;
EV_x = V'*WVx*V;
EA_x = A'*WAx*A;
E_x = EX_x + EV_x + EA_x + T;
G2X_x = X'*WXx*x;
G2V_x = V'*WVx*x_t;
G2A_x = A'*WAx*x_tt;
G2_x = G2X_x + G2V_x + G2A_x;


EX_y = X'*WXy*X;
EV_y = V'*WVy*V;
EA_y = A'*WAy*A;
E_y = EX_y + EV_y + EA_y + T;
G2X_y = X'*WXy*y;
G2V_y = V'*WVy*y_t;
G2A_y = A'*WAy*y_tt;
G2_y = G2X_y + G2V_y + G2A_y;

withConstraints = 0;

if withConstraints == 0
    m_x = E_x\G2_x;
    m_y = E_y\G2_y;

    Cm_x = inv(E_x);
    Cm_y = inv(E_y);
else
    C_x = zeros(M,M); % MxM
    C_y = zeros(M,M); % MxM

    G1 = [E_x,C_x,F',zeros(M,NC);
          C_y,E_y,zeros(M,NC),F';
          F,zeros(NC,M),zeros(NC,NC),zeros(NC,NC);
          zeros(NC,M),F,zeros(NC,NC),zeros(NC,NC)];

    G2 = [G2_x;G2_y;hx';hy'];
    m = G1\G2;
    m_x = m(1:M);
    m_y = m(M+1:2*M);

    Cm_x = inv(E_x) - inv(E_x)*F'*inv(F*inv(E_x)*F')*F*inv(E_x);
    Cm_y = inv(E_y) - inv(E_y)*F'*inv(F*inv(E_y)*F')*F*inv(E_y);
end

end

function [m_x,m_y,Cm_x,Cm_y] = ComputeSolutionWithElim( X, V, A, T, WXx, WXy, Sigma_x, Sigma_y, Sigma_xx, Sigma_yy, F, hx, hy, f0, M, NC, x, y, x_t,y_t, x_tt, y_tt )

includeX = 1;
includeV = 0;
includeA = 0;

EX_x = X'*WXx*X;
EV_x = V'*(Sigma_x\V);
EA_x = A'*(Sigma_xx\A);
E_x = includeX*EX_x + includeV*EV_x + includeA*EA_x + T;
G2X_x = X'*WXx*x;
G2V_x = V'*(Sigma_x\x_t);
G2A_x = A'*(Sigma_xx\x_tt);
G2_x = includeX*G2X_x + includeV*G2V_x + includeA*G2A_x;


EX_y = X'*WXy*X;
EV_y = V'*(Sigma_y\V);
EA_y = A'*(Sigma_yy\A);
E_y = includeX*EX_y + includeV*EV_y + includeA*EA_y + T;
G2X_y = X'*WXy*y;
G2V_y = V'*(Sigma_y\y_t);
G2A_y = A'*(Sigma_yy\y_tt);
G2_y = includeX*G2X_y + includeV*G2V_y + includeA*G2A_y;

withConstraints = 0;

if withConstraints == 0
    m_x = E_x\G2_x;
    m_y = E_y\G2_y;

    Cm_x = inv(E_x);
    Cm_y = inv(E_y);
else
    C_x = zeros(M,M); % MxM
    C_y = zeros(M,M); % MxM

    G1 = [E_x,C_x,F',zeros(M,NC);
          C_y,E_y,zeros(M,NC),F';
          F,zeros(NC,M),zeros(NC,NC),zeros(NC,NC);
          zeros(NC,M),F,zeros(NC,NC),zeros(NC,NC)];

    G2 = [G2_x;G2_y;hx';hy'];
    m = G1\G2;
    m_x = m(1:M);
    m_y = m(M+1:2*M);

    Cm_x = inv(E_x) - inv(E_x)*F'*inv(F*inv(E_x)*F')*F*inv(E_x);
    Cm_y = inv(E_y) - inv(E_y)*F'*inv(F*inv(E_y)*F')*F*inv(E_y);
end

end


