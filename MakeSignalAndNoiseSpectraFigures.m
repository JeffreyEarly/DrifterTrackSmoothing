addpath('./support');
% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

load('sample_data/SyntheticTrajectories.mat')

timescale = 60;
indices = find( t/timescale >= 0 & t/timescale <= 42);

markersize  = 2*scaleFactor;
range = indices(1:10:end);

FigureSize = [50 50 figure_width_large+7 300*scaleFactor];

fig1 = figure('Units', 'points', 'Position', FigureSize)
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];



subplot(2,1,1)
errorbar(t(range)/timescale,y(range),position_error*ones(size(t(range))),'ko', 'LineWidth', 0.5*scaleFactor, 'MarkerSize', markersize^2, 'MarkerFaceColor', 'w')
hold on
xlim([min(t(range))/timescale max(t(range))/timescale]);
ylim([1.2*min(y(range)) 1.2*max(y(range))])

range = indices(1:100:end);
errorbar(t(range)/timescale,y(range),position_error*ones(size(t(range))),'ko', 'LineWidth', 1.5*scaleFactor, 'MarkerSize', markersize^2, 'MarkerFaceColor', 'k')

xlabel('minutes', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
ylabel('meters', 'FontSize', figure_axis_label_size, 'FontName', figure_font);



D = FiniteDifferenceMatrixNoBoundary(1,t,1);



dt = t(2)-t(1);
cv = D*(x + sqrt(-1)*y);
cepsilon = D*(epsilon_x + sqrt(-1)*epsilon_y);
[psi,lambda]=sleptap(size(cv,1));
[f,spp,snn,spn]=mspec(dt,cv,psi);
[f,spp_e,snn_e,spn_e]=mspec(dt,cepsilon,psi);

% sqrt(2)*position_error/dt
% std(imag(cepsilon))

ylimit = [1e-4 4e2];

subplot(2,1,2)
plot(f*timescale,vmean([snn, spp],2), 'k', 'LineWidth', 2)
hold on
plot(f*timescale,vmean([snn_e, spp_e],2), 'r', 'LineWidth', 2)
xlog, ylog
xlim([min(f*timescale) max(f*timescale)])
ylim(ylimit)

dt_gamma1 = 3*position_error/sigma_u;
f_gamma1 = timescale/dt_gamma1;
% plot([f_gamma1 f_gamma1],ylimit, 'LineWidth', 0.5*scaleFactor, 'Color', 0.4*[1.0 1.0 1.0]);
% 
dt_gamma10 = 3*position_error/(sigma_u*10);
f_gamma10 = timescale/dt_gamma10;
% plot([f_gamma10 f_gamma10],ylimit, 'LineWidth', 0.5*scaleFactor, 'Color', 0.4*[1.0 1.0 1.0]);
% 
dt_gamma01 = 3*position_error/(sigma_u*0.1);
f_gamma01 = timescale/dt_gamma01;
% plot([f_gamma01 f_gamma01],ylimit, 'LineWidth', 0.5*scaleFactor, 'Color', 0.4*[1.0 1.0 1.0]);

xlabel('cycles per minute', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
ylabel('power (m^2/s)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

print('-depsc2', 'figures/synthetic_process_and_spectrum.eps')

dof1 = 1 + 3*position_error/(sigma_u*dt);
dof10 = 1 + 3*position_error/(sigma_u*10*dt);
dof100 = 1 + 3*position_error/(sigma_u*100*dt);

return

indices = 1:1:floor(length(t)/5);
indicesAll = 1:max(indices);
x_obs = x(indices) + epsilon_x(indices);
y_obs = y(indices) + epsilon_y(indices);
t_obs = t(indices);
sigma = position_error;
S = 3;
T = 2;
[m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t_obs,x_obs,y_obs,sigma,S,T);

D_obs = FiniteDifferenceMatrixNoBoundary(1,t_obs,1);
X = squeeze(B(:,:,1));
cv_tension = D_obs*(X*m_x + sqrt(-1)*X*m_y);
cepsilon_tension = D_obs*( (x_obs-X*m_x) + sqrt(-1)*(y_obs-X*m_y));
[psi,lambda]=sleptap(size(cv_tension,1));
[f,spp,snn,spn]=mspec(t_obs(2)-t_obs(1),cv_tension,psi);
[f,spp_e,snn_e,spn_e]=mspec(dt,cepsilon_tension,psi);
plot(f*timescale,vmean([snn, spp],2), 'g', 'LineWidth', 2)
plot(f*timescale,vmean([snn_e, spp_e],2), 'b', 'LineWidth', 2)

fig1 = figure, plot(t_obs,X*m_x, 'g'), hold on, plot(t,x,'k')

% packfig(maxK,maxD)
fig1 = tightfig;
fig1.Position = FigureSize;
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];
fig1.PaperPositionMode = 'auto';



a_vals = 10.^(linspace(-2,1,15))';
rms_error = zeros(size(a_vals));
for i = 1:length(a_vals)
   [m_x,m_y,Cm_x,Cm_y,B,Bq,tq] = smooth_interpolate_gaussian_noise(t_obs,x_obs,y_obs,sigma,S,T,a_vals(i));
   X = squeeze(B(:,:,1));
   rms_error(i) = std( (x(indices)-X*m_x) + sqrt(-1)*(y(indices)-X*m_y) );   
end

figure
scatter(a_vals, rms_error);
