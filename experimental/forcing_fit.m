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
function [m_x,m_y,Cm_x,Cm_y,X,V,A,J] = forcing_fit(t,x,y,dx,dy,T_decorrelation,M,S,v0,lat0, weight_function)

if (length(t) ~= length(x) || length(t) ~= length(y) )
   disp('The time series are not consistent lengths');
   return;
end

if S == 3
    addpath('./cubic_splines');
    spline_int = @(t) cspline_int(t);
    spline = @(t) cspline(t);
    spline_t = @(t) cspline_t(t);
    spline_tt = @(t) cspline_tt(t);
    spline_ttt = @(t) cspline_ttt(t);
elseif S == 5
    addpath('./quintic_splines');
    spline_int = @(t) qspline_int(t);
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
Xi = zeros(N,M);
X = zeros(N,M);
V = zeros(N,M);
A = zeros(N,M);
J = zeros(N,M);
for i=1:N
    for j=1:M
        t_norm=(t(i)-t(1))/t_knot - (j - 1 - floor(S/2));
        Xi(i,j)=spline_int(t_norm);
        X(i,j)=spline(t_norm);
        V(i,j)=spline_t(t_norm);
        A(i,j)=spline_tt(t_norm);
        J(i,j)=spline_ttt(t_norm);
    end
end
Xi = Xi*t_knot;
V = V/t_knot;
A = A/(t_knot^2);
J = J/(t_knot^3);

% set up D matrix
Xq = zeros(Q,M); % Same as the A matrix above, but on the quadrature (q) grid.
Vq = zeros(Q,M);
for q=1:Q
    tq = t(1) + (q-1)*DT;
    for j=1:M
        t_norm=(tq-t(1))/t_knot - (j - 1 - floor(S/2));
        Xq(q,j)=spline(t_norm);
        Vq(q,j)=spline_t(t_norm);
    end
end
Vq = Vq/t_knot;

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

% 2nd order 2nd derivative matrix, no boundary terms
Diff2 = FiniteDifferenceMatrix(2,t,2,2,2);
Diff2 = Diff2(2:(N-1),:);

Sigma_x=zeros(N-2,N-2);
Sigma_y=zeros(N-2,N-2);
for i=1:size(Sigma_x,1)
    for j=1:size(Sigma_x,2)
        Sigma_x(i,j) = sum(Diff2(i,:).*Diff2(j,:).*dx'.*dx');
        Sigma_y(i,j) = sum(Diff2(i,:).*Diff2(j,:).*dy'.*dy');
    end
end

WAx = inv(Sigma_x);
WAy = inv(Sigma_y);

WXx = diag(1./(dx.^2));
WXy = diag(1./(dy.^2));

% Acceleration, but w/o the first and last point
Af = A(2:(N-1),:);
% Af = A;

T = gamma*(Vq'*Vq);

EX_x = X'*WXx*X;
EA_x = Af'*WAx*Af;
E_x = EX_x + EA_x + T;
G2X_x = X'*WXx*x;
G2A_x = Af'*WAx*Diff2*x;
G2_x = G2X_x + G2A_x;


EX_y = X'*WXy*X;
EA_y = Af'*WAy*Af;
E_y = EX_y + EA_y + T;
G2X_y = X'*WXy*y;
G2A_y = Af'*WAy*Diff2*y;
G2_y = G2X_y + G2A_y;

withConstraints = 1;

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


