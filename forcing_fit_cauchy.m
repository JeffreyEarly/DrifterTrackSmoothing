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
function [m_x,m_y,Cm_x,Cm_y,X,V,A,J] = forcing_fit_cauchy(t,x,y,dx,dy,T_decorrelation,M,S,v0,lat0, weight_function)

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

% Rows are the N observations
% Columns are the M splines
t_norm = zeros(N,M);
for i=1:N
    for j=1:M
        t_norm(i,j)=(t(i)-t(1))/t_knot - (j - 1 - floor(S/2));   
    end
end
X=spline(t_norm);
V=spline_t(t_norm)/t_knot;
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
Diff2 = FiniteDifferenceMatrix(2,t,2,2,10);
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

% Tension
T = gamma*(Vq'*Vq);
x_t = Diff2*x;
y_t = Diff2*y;

[m_x,m_y,Cm_x,Cm_y] = ComputeSolution( X, V, Af, T, WXx, WXy, WAx, WAy, F, hx, hy, f0, M, NC, x, y, x_t, y_t );

error_y_previous = dy;
rel_error = 1.0;
repeats = 1;
while (rel_error > 0.01)
    dx1 = X*m_x - x;
    dy1 = X*m_y - y;
    
    dx1 = max(0.00001*ones(size(dx1)),abs(dx1));
    dy1 = max(0.00001*ones(size(dy1)),abs(dy1));
%     ds1 = max(abs(dx1),abs(dy1));
    
    dx2 = dx1./weight_function(dx1./dx);
    dy2 = dy1./weight_function(dy1./dx);
%     ds2 = ds1./weight_function(ds1./dx);
    
%     dx2 = ds2;
%     dy2 = ds2;
    
    Sigma_x=zeros(N-2,N-2);
    Sigma_y=zeros(N-2,N-2);
    for i=1:size(Sigma_x,1)
        for j=1:size(Sigma_x,2)
            Sigma_x(i,j) = sum(Diff2(i,:).*Diff2(j,:).*dx2'.*dx2');
            Sigma_y(i,j) = sum(Diff2(i,:).*Diff2(j,:).*dy2'.*dy2');
        end
    end

%     WAx = inv(Sigma_x);
%     WAy = inv(Sigma_y);

    WXx = diag(1./(dx2.^2));
    WXy = diag(1./(dy2.^2));
    
    %         dbstop if warning
%     [m_x,m_y,Cm_x,Cm_y] = ComputeSolution( X, V, Af, T, WXx, WXy, WAx, WAy, F, hx, hy, f0, M, NC, x, y, x_t, y_t );
    [m_x,m_y,Cm_x,Cm_y] = ComputeSolutionWithElim( X, V, Af, T, WXx, WXy, Sigma_x, Sigma_y, F, hx, hy, f0, M, NC, x, y, x_t, y_t );
    
    rel_error = max((dx2-error_y_previous)./dx2);
    error_y_previous=dx2;
    repeats = repeats+1;
    
    if (repeats == 100)
        disp('Failed to converge after 100 iterations.');
        break;
    end
    
end

end


function [m_x,m_y,Cm_x,Cm_y] = ComputeSolution( X, V, A, T, WXx, WXy, WAx, WAy, F, hx, hy, f0, M, NC, x, y, x_t, y_t )

EX_x = X'*WXx*X;
EA_x = A'*WAx*A;
E_x = EX_x + EA_x + T;
G2X_x = X'*WXx*x;
G2A_x = A'*WAx*x_t;
G2_x = G2X_x + G2A_x;


EX_y = X'*WXy*X;
EA_y = A'*WAy*A;
E_y = EX_y + EA_y + T;
G2X_y = X'*WXy*y;
G2A_y = A'*WAy*y_t;
G2_y = G2X_y + G2A_y;

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

function [m_x,m_y,Cm_x,Cm_y] = ComputeSolutionWithElim( X, V, A, T, WXx, WXy, Sigma_x, Sigma_y, F, hx, hy, f0, M, NC, x, y, x_t, y_t )

EX_x = X'*WXx*X;
EA_x = A'*(Sigma_x\A);
E_x = EX_x + EA_x + T;
G2X_x = X'*WXx*x;
G2A_x = A'*(Sigma_x\x_t);
G2_x = G2X_x + G2A_x;


EX_y = X'*WXy*X;
EA_y = A'*(Sigma_y\A);
E_y = EX_y + EA_y + T;
G2X_y = X'*WXy*y;
G2A_y = A'*(Sigma_y\y_t);
G2_y = G2X_y + G2A_y;

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


