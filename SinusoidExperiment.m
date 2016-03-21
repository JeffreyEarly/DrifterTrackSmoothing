% We are going to create a Square Wave, with the number of points at each
% crest and trough specified by DOF.

DOF = 5;

% Create a signal
dt = 0.5/DOF;
t=(0:dt:20)';


x=cos(2*pi*t);
x = 2*(floor(0.5*cos(2*pi*t+pi/4))+0.5);
% x = 2*(round(0.5*sin(2*pi*t)+0.5)-0.5);
u = diff(x)/dt;
a = diff(u)/dt;

% Compute its spectrum
[f,s]=mspec(dt,x,[]);
[f_u,s_u]=mspec(dt,u,[]);
[f_a,s_a]=mspec(dt,a,[]);

df = f(2)-f(1);

sigma = 1.0;
epsilon = sigma*randn(size(t));
epsilon_t = diff(epsilon)/dt;
epsilon_tt = diff(epsilon_t)/dt;

[~,s_e]=mspec(dt,epsilon,[]);
[~,s_e_u]=mspec(dt,epsilon_t,[]);
[~,s_e_a]=mspec(dt,epsilon_tt,[]);


% Check that the variances match (ignoring the mean for the moment)
x_variance = std(x).^2;
s_variance = sum(s)*f(2)/(2*pi);
fprintf('x variance: %g, spectrum variance: %g\n', x_variance, s_variance);
u_variance = std(u).^2;
s_u_variance = sum(s_u)*f_u(2)/(2*pi);
fprintf('u variance: %g, spectrum variance: %g\n', u_variance, s_u_variance);
a_variance = std(a).^2;
s_a_variance = sum(s_a)*f_a(2)/(2*pi);
fprintf('a variance: %g, spectrum variance: %g\n', a_variance, s_a_variance);

figure
subplot(2,1,1)
plot(t,x)
ylim([-1.1 1.1])
subplot(2,1,2)
plot(f/(2*pi),s)
hold on
plot(f_u/(2*pi),s_u)
plot(f_a/(2*pi),s_a)
plot(f/(2*pi),s_e,'k')
plot(f/(2*pi),sigma*sigma*dt*ones(size(f)),'k')
plot(f/(2*pi),s_e_u,'k')
plot(f/(2*pi),sigma*sigma*dt*f.*f,'k')
legend('x', 'u', 'a')

plot(f/(2*pi),s.*f.*f, 'r')

xlim([min(f) max(f)]/(2*pi))
ylim([1e-5 10*max(s_a)])
xlog,ylog