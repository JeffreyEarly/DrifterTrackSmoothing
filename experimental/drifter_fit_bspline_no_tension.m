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
function [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = drifter_fit_bspline_no_tension(t,x,y,dx,dy,S,t_knot,weight_function)

if (length(t) ~= length(x) || length(t) ~= length(y) )
   disp('The time series are not consistent lengths');
   return;
end

B = bspline(t,t_knot,S+1);
X = squeeze(B(:,:,1));
N = length(y); % also size(X,1);
M = size(X,2); % number of splines

% Now we need a quadrature (integration) grid that is finer
Q = 10*N; % number of points on the quadrature grid
tq = linspace(t(1),t(end),Q)';
Bq = bspline(tq,t_knot,S+1);

% set up F matrix and h vector for constraints
NC = 2;
F=zeros(NC,M);
F(1,:) = squeeze(B(1,:,end));
F(2,:) = squeeze(B(end,:,end));
h(1) = 0.;
h(2) = 0.;

Wx = diag(1./(dx.^2));
Wy = diag(1./(dy.^2));

dbstop if warning
[m_x,m_y,Cm_x,Cm_y] = ComputeSolution( X, Bq, F, F, Wx, Wy, M, NC, x, y, h );

error_x_previous = dx;
error_y_previous = dy;
rel_error = 1.0;
repeats = 1;
while (rel_error > 0.01)
    dx2 = weight_function(X*m_x - x);
    dy2 = weight_function(X*m_y - y);
    
    Wx = diag(1./(dx2));
    Wy = diag(1./(dy2));
    
    [m_x,m_y,Cm_x,Cm_y] = ComputeSolution( X, Bq, F, F, Wx, Wy, M, NC, x, y, h );
    
    rel_error = max(max( (dx2-error_x_previous)./dx2 ),max( (dy2-error_y_previous)./dy2 ));
    error_x_previous=dx2;
    error_y_previous=dy2;
    repeats = repeats+1;
    
    if (repeats == 100)
        disp('Failed to converge after 100 iterations.');
        break;
    end
end

dx=X*m_x - x;
dy=X*m_y - y;
a=dx'*Wx*dx/N;
b=dy'*Wy*dy/N;
chi = 2/(a+b);
sigma2 = sum(dx.*dx + dy.*dy)/(2*N);
w_rms = median(sqrt(1./[diag(Wx);diag(Wy)]));
fprintf('sum(dx/sigma)^2=%f, sum(dy/sigma)^2=%f, chi=%f, sigma = %f, w_rms=%f\n',a,b,chi,sqrt(sigma2),w_rms);


end


function [m_x,m_y,Cm_x,Cm_y] = ComputeSolution( X, Bq, F, Fy, Wx, Wy, M, NC, x, y, h )
% A (and D) matrix:
% Rows are the N observations
% Columns are the M splines
% F is NCxM
% H is NCx1


f0=0;
% set up inverse theory matrices
E_x = X'*Wx*X; % MxM
E_y = X'*Wy*X; % MxM

    m_x = E_x\(X'*Wx*x);
    m_y = E_y\(X'*Wy*y);

    Cm_x = inv(E_x);
    Cm_y = inv(E_y);
return
Xq = squeeze(Bq(:,:,1));
Vq = squeeze(Bq(:,:,2));
C_x = -f0*0*(Vq'*Xq); % MxM
C_y = -f0*0*(Xq'*Vq); % MxM

    G1 = [E_x,C_x,F',zeros(M,NC);
          C_y,E_y,zeros(M,NC),F';
          F,zeros(NC,M),zeros(NC,NC),zeros(NC,NC);
          zeros(NC,M),F,zeros(NC,NC),zeros(NC,NC)];
      
G2 = [X'*Wx*x;X'*Wy*y;h';h'];
m = G1\G2;
m_x = m(1:M);
m_y = m(M+1:2*M);

% model parameter covariance matrix
%Cm_x = inv(E_x) - inv(E_x)*F'*inv(F*inv(E_x)*F')*F*inv(E_x);
%Cm_y = inv(E_y) - inv(E_y)*F'*inv(F*inv(E_y)*F')*F*inv(E_y);
Cm_x = zeros(size(E_x));
Cm_y = zeros(size(E_x));
end

