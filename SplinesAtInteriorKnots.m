K = 4;
S = K-1;
t_knot = linspace(0,10,11)';
t_knot(3) = [];
t_knot(5) = [];

t = linspace(-1,11,121)';
% t = 10;

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
% pre_knots = t_knot(1) + (t_knot(2)-t_knot(1))*(-floor(S/2):-1)';
% post_knots = t_knot(end) + (t_knot(end)-t_knot(end-1))*(1:floor(S/2))';
%t_knot = [pre_knots; t_knot; post_knots];
t_knot = [repmat(t_knot(1),S,1); t_knot; repmat(t_knot(end),S,1)];
dt_knot = diff(t_knot);
t_knot2 = t_knot + [dt_knot; dt_knot(end)];

% numer of knots
M = length(t_knot);

% This is true assuming the original t_knot was strictly monotonically
% increasing (no repeat knots) and we added repeat knots at the beginning
% and end of the sequences.
N_splines = M - 2*S + 2*floor(S/2);

% number of collocation points
N = length(t);

% Rows are the N collocation points
% Columns are the M splines
X = zeros(N,N_splines,K); % This will contain all splines and their derivatives
XB = zeros(N,N_splines,K); % This will contain all splines through order K
for t_i=1:N % loop through all N collocation points
    i = find( t_knot <= t(t_i) & t(t_i) < t_knot2, 1, 'last' );
    if isempty(i)
        if t(t_i) < t_knot(1)
            i = find( t_knot == t_knot(1), 1, 'last');
            continue; %This continue means we don't need to set b(1) = 0; or check indices on the delta_r line
        elseif t(t_i) == t_knot(end)
            i = find( t_knot < t(t_i), 1, 'last'); 
        else
            i = find( t_knot < t_knot(end), 1, 'last');
            continue; %b(1) = 0;
        end
    end
    
    delta_r = zeros(K,1);
    delta_l = zeros(K,1);
    
    XB(t_i,i,1) = 1;
    
    b = zeros(K,1); b(1) = 1;
    for j=1:(K-1) % loop through splines of increasing order
       delta_r(j) = t_knot(i+j) - t(t_i);
       delta_l(j) = t(t_i) - t_knot(i+1-j);
       
       saved = 0;
       for r=1:j % loop through the nonzero splines
           term = b(r)/(delta_r(r) + delta_l(j+1-r));
           b(r) = saved + delta_r(r)*term;
           saved = delta_l(j+1-r)*term;
       end
       b(j+1) = saved;
       
       indices = max(1,i-j):i;
       XB(t_i,indices,j+1) = b(1:length(indices));
    end
    
    indices = max(1,i-K+1):i;
    X(t_i,indices,1) = b(1:length(indices));
    
end

% diff_coeff = @(a,r,m) (K-m)*(a(r)-a(r-1))/(t_knot(r+K-m) - t_knot(r));
diff_coeff = @(a,r,m) (K-m)*(a(2)-a(1))/(t_knot(r+K-m) - t_knot(r));

for r=1:N_splines
    % alpha mimics equation X.16 in deBoor's PGS, but localized to avoid
    % the zero elements.
    alpha = zeros(S+2,S+2); % row is the coefficient, column is the derivative (1=0 derivatives)
    alpha(2,1) = 1;
    for m=1:S
        for i=1:(m+1)
            a = alpha(:,m);
            alpha(i+1,m+1) = diff_coeff(a(i:end),r+i-1,m);
            if isinf(alpha(i+1,m+1)) || isnan(alpha(i+1,m+1))
                alpha(i+1,m+1) = 0;
            end
            if r+i-1>N_splines
                B = zeros(N,1);
            else
                B = XB(:,r+i-1,K-m);
            end
            X(:,r,m+1) = X(:,r,m+1) + alpha(i+1,m+1)*B;
        end
    end
end

% for i=1:N_splines
%     alpha = zeros(m+1,m);
%     for m=1:S
%         delta_l = t_knot(i+K-m) - t_knot(i);
%         delta_r = t_knot(i+K-m+1) - t_knot(i+1);
%         if i+1>N_splines
%             Br = zeros(N,1);
%         else
%             Br = XB(:,i+1,K-m);
%         end
%         X(:,i,m+1) = (K-m)*(XB(:,i,K-m)/delta_l - Br/delta_r);
%     end
% end

k=1;
figure, plot(t,X(:,k,1),'LineWidth', 2), hold on
plot(t,X(:,k,2),'LineWidth', 2)
plot(t,X(:,k,3),'LineWidth', 2)
plot(t,X(:,k,4),'LineWidth', 2)
% plot(t,DB(:,k,K-3))

plot(t,vdiff(t(2)-t(1),X(:,k,1),1))
plot(t,vdiff(t(2)-t(1),vdiff(t(2)-t(1),X(:,k),1),1))
% ylim([min(X(:,k))*1.05 max(X(:,k))*1.05])
legend('X', 'X_t', 'X_{tt}', 'X_{ttt}', 'diff(X)', 'diff(diff(X))')
vlines(t_knot,'g--')

figure
plot(t,X(:,:,1))
ylim([min(min(X(:,:,1)))*1.05 max(max(X(:,:,1)))*1.05])
vlines(t_knot,'g--')