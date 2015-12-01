S = 3;
t_knot = linspace(0,10,11)';
t_knot(3) = [];
t_knot(5) = [];

t = linspace(0,10,101)';

% if S == 3
%     addpath('./cubic_splines');
%     spline = @(t) cspline(t);
%     spline_t = @(t) cspline_t(t);
%     spline_tt = @(t) cspline_tt(t);
%     spline_ttt = @(t) cspline_ttt(t);
% elseif S == 5
%     addpath('./quintic_splines');
%     spline = @(t) qspline(t);
%     spline_t = @(t) qspline_t(t);
%     spline_tt = @(t) qspline_tt(t);
%     spline_ttt = @(t) qspline_ttt(t);
% else
%     disp('Whoops! I only know how to deal with splines of order 3 and 5.')
%     return;
% end

% Now we need to add knot points past the beginning and end.
%
% M is given as the number of splines on the interval [t(1), t(end)]. We
% need to add extra splines before an after this, depending on order of the
% spline.
pre_knots = t_knot(1) + (t_knot(2)-t_knot(1))*(-floor(S/2):-1)';
post_knots = t_knot(end) + (t_knot(end)-t_knot(end-1))*(1:floor(S/2))';
t_knot = [pre_knots; t_knot; post_knots];
dt_knot = diff(t_knot);

% numer of knots
M = length(t_knot);

% number of collocation points
N = length(t);

% Rows are the N collocation points
% Columns are the M splines
X = zeros(N,M);
for t_i=1:N % loop through all N collocation points
    i = find(t(t_i)<t_knot,1,'first')-1;
    if isempty(i)
        i = M;
    end
    
    delta_r = zeros(S,1);
    delta_l = zeros(S,1);
    b = zeros(S,1); b(1) = 1;
    
    for j=1:(S-1) % loop through splines of increasing order
        if (i+1-j) < 1 || i+j > M
            continue;
        end
       delta_r(j) = t_knot(i+j) - t(t_i);
       delta_l(j) = t(t_i) - t_knot(i+1-j);
       
       saved = 0;
       for r=1:j % loop through the nonzero splines
           term = b(r)/(delta_r(r) + delta_l(j+1-r));
           b(r) = saved + delta_r(r)*term;
           saved = delta_l(j+1-r)*term;
       end
       b(j+1) = saved;
    end
    
    indices = max(1,i-S+1):i;
    X(t_i,indices) = b(1:length(indices));
    
end

k=1;
figure, plot(t,X(:,k)),ylim([min(X(:,k))*1.05 max(X(:,k))*1.05]),vlines(t_knot,'g--')

figure
plot(t,X)
ylim([min(min(X))*1.05 max(max(X))*1.05])
vlines(t_knot,'g--')