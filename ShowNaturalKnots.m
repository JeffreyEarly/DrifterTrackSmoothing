t = (0:6)';
x = sin(2*pi*t/4);
sigma = 1;
dx = sigma*ones(size(x));

N = length(t);

K = 1;
S = K-1;
[t_knot] = NaturalKnotsForSpline( t, K );
w = @(z)(sigma*sigma);
[m_x,Cm_x,B] = bspline_fit_no_tension(t,x,dx,S,t_knot,w);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Position Figures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

tq = linspace(t(1),t(end),10*N)';
Bq = bspline(tq,t_knot,S+1);

subplot(3,1,1)
plot(tq,squeeze(Bq(:,:,1))*m_x,'b','LineWidth',1.5), hold on
scatter(t,squeeze(B(:,:,1))*m_x,6^2,'filled')
vlines(t_knot0, 'g--')

return;

subplot(3,1,2)
plot(t,x_true,'k','LineWidth',1), hold on
plot(tq1,x_fit1,'b','LineWidth',1.5)
scatter(t,x,6^2,'filled')
vlines(t_knot1, 'g--')

subplot(3,1,3)
plot(t,x_true,'k','LineWidth',1), hold on
plot(tq2,x_fit2,'b','LineWidth',1.5)
scatter(t,x,6^2,'filled')
vlines(t_knot2, 'g--')
