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
function [m_x,m_y,Cm_x,Cm_y,X,V,A,J,Xq,Vq] = drifter_fit_bspline(t,x,y,dx,dy,S,v0, weight_function)

if (length(t) ~= length(x) || length(t) ~= length(y) )
   disp('The time series are not consistent lengths');
   return;
end

B = bspline(t,t,S+1);
X = squeeze(B(:,:,1));
V = squeeze(B(:,:,2));
A = squeeze(B(:,:,3));
J = squeeze(B(:,:,4));
N = length(y); % also size(X,1);
M = size(X,2); % number of splines

% Now we need a quadrature (integration) grid that is finer
Q = 10*N; % number of points on the quadrature grid
tq = linspace(t(1),t(end),Q)';
Bq = bspline(tq,t,S+1);
Xq = squeeze(Bq(:,:,1));
Vq = squeeze(Bq(:,:,2));

gamma = 1/(v0*v0*Q);

% set up F matrix and h vector for constraints
NC = 2;
F=zeros(NC,M);
F(1,:) = A(1,:);
F(2,:) = A(end,:);
h(1) = 0.;
h(2) = 0.;

Wx = diag(1./(dx.^2));
Wy = diag(1./(dy.^2));

dbstop if warning
[m_x,m_y,Cm_x,Cm_y] = ComputeSolution( X, Xq, Vq, F, F, Wx, Wy, gamma, M, NC, x, y, h );

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
    
    Wx = diag(1./(dx2.^2));
    Wy = diag(1./(dy2.^2));
    Ws = diag(1./(ds2.^2));
    
    %          dbstop if warning
    [m_x,m_y,Cm_x,Cm_y] = ComputeSolution( X, Xq, Vq, F, F, Wx, Wy, gamma, M, NC, x, y, h );
    
    rel_error = max((dx2-error_y_previous)./dx2);
    error_y_previous=dx2;
    repeats = repeats+1;
    
    if (repeats == 100)
        disp('Failed to converge after 100 iterations.');
        break;
    end
end




end


function [m_x,m_y,Cm_x,Cm_y] = ComputeSolution( X, Xq, Vq, F, Fy, Wx, Wy, gamma, M, NC, x, y, h )
% A (and D) matrix:
% Rows are the N observations
% Columns are the M splines
% F is NCxM
% H is NCx1
f0=0;
% set up inverse theory matrices
E_x = X'*Wx*X + gamma*(Vq'*Vq); % MxM
E_y = X'*Wy*X + gamma*(Vq'*Vq); % MxM
C_x = -f0*gamma*(Vq'*Xq); % MxM
C_y = -f0*gamma*(Xq'*Vq); % MxM

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

