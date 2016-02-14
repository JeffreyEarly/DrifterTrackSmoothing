
TotalPeriods = 10;
Amplitude = 4;
T = 100;
sigma = 1;
t=linspace(0,TotalPeriods*T,200)';

x_analytic = @(t) Amplitude*sin(2*pi*t/T);
v_analytic = @(t) Amplitude*(2*pi/T)*cos(2*pi*t/T);
a_analytic = @(t) -Amplitude*(2*pi/T)^2*sin(2*pi*t/T);

% create the 'true' signal
x_true = x_analytic(t);

% Create some Gaussian noise
epsilon = sigma*randn(size(x_true));

% Contaminate the signal
x = x_true + epsilon;

% Define the noise pdfs
p = @(z) exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi));
w = @(z)(sigma*sigma);

T = 2; a = 1e-2;

S = T+1; K = S+1;


tension = zeros(S,1);
tension(T) = 1/a^2;
[m_x,Cm_x,B,Bq,tq] = bspline_fit_with_tension(t,x,ones(size(x))*sigma,S,tension,w);
x_error = x - squeeze(B(:,:,1))*m_x;
x_fit = squeeze(Bq(:,:,1))*m_x;
v_fit = squeeze(Bq(:,:,2))*m_x;

std(x_error)

figure
plot(t,x_true), hold on
plot(tq,x_fit,'g')
scatter(t,x,5)

figure
histogram(v_analytic(t)), hold on
histogram(v_fit)