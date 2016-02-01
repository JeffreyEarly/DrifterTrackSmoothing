t = (0:10)';
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

totalPlots = 5;
subplot(totalPlots,1,1)
plot(tq,squeeze(Bq(:,:,1))*m_x,'b','LineWidth',1.5), hold on
scatter(t,squeeze(B(:,:,1))*m_x,6^2,'filled')
vlines(t_knot, 'g--')


K = 2;
S = K-1;
[t_knot] = NaturalKnotsForSpline( t, K );
[m_x,Cm_x,B] = bspline_fit_no_tension(t,x,dx,S,t_knot,w);
Bq = bspline(tq,t_knot,S+1);


subplot(totalPlots,1,2)
plot(tq,squeeze(Bq(:,:,1))*m_x,'b','LineWidth',1.5), hold on
scatter(t,squeeze(B(:,:,1))*m_x,6^2,'filled')
vlines(t_knot, 'g--')




K = 3;
S = K-1;
[t_knot] = NaturalKnotsForSpline( t, K );
[m_x,Cm_x,B] = bspline_fit_no_tension(t,x,dx,S,t_knot,w);
Bq = bspline(tq,t_knot,S+1);

subplot(totalPlots,1,3)
plot(tq,squeeze(Bq(:,:,1))*m_x,'b','LineWidth',1.5), hold on
scatter(t,squeeze(B(:,:,1))*m_x,6^2,'filled')
vlines(t_knot, 'g--')

K = 4;
S = K-1;
[t_knot] = NaturalKnotsForSpline( t, K );
[m_x,Cm_x,B] = bspline_fit_no_tension(t,x,dx,S,t_knot,w);
Bq = bspline(tq,t_knot,S+1);

subplot(totalPlots,1,4)
plot(tq,squeeze(Bq(:,:,1))*m_x,'b','LineWidth',1.5), hold on
scatter(t,squeeze(B(:,:,1))*m_x,6^2,'filled')
vlines(t_knot, 'g--')

K = 5;
S = K-1;
[t_knot] = NaturalKnotsForSpline( t, K );
[m_x,Cm_x,B] = bspline_fit_no_tension(t,x,dx,S,t_knot,w);
Bq = bspline(tq,t_knot,S+1);

subplot(totalPlots,1,5)
plot(tq,squeeze(Bq(:,:,1))*m_x,'b','LineWidth',1.5), hold on
scatter(t,squeeze(B(:,:,1))*m_x,6^2,'filled')
vlines(t_knot, 'g--')