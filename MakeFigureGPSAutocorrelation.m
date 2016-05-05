% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults
addpath('support')

stride = 1;
maxT = 1.5*3600;

load('sample_data/motionless_garmin_epix.mat')
x=x-mean(x);
y=y-mean(y);

x_out = x(1:stride:end);
y_out = y(1:stride:end);

dt = (0:stride:maxT)';
n = (0:(length(dt)-1))';
[ACx, DOFx] = Autocorrelation(x_out, length(dt)-1);
[ACy, DOFy] = Autocorrelation(y_out, length(dt)-1);

AC = ACx + ACy;
DOF = DOFx + DOFy;
 
load('sample_data/motionless_garmin_edge_705.mat')

x=x-mean(x);
y=y-mean(y);

x_out = x(1:stride:end);
y_out = y(1:stride:end);

dt = (0:stride:maxT)';
n = (0:(length(dt)-1))';
[ACx, DOFx] = Autocorrelation(x_out, length(dt)-1);
[ACy, DOFy] = Autocorrelation(y_out, length(dt)-1);

AC = AC + ACx + ACy;
DOF = DOF + DOFx + DOFy;

AC = AC/4;

SE_indep = dt(2:end);
SE =  sqrt((1 + 2*cumsum(AC.^2))./DOF);
SE(1) = sqrt(1/DOF(1)); % first point is lag 1
SE(end) = []; % there is no end point

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Autocorrelation figures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FigureSize = [50 50 figure_width_1col+4 150*scaleFactor];

fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

s = 1/60;

plot(s*dt,AC, 'LineWidth',1*scaleFactor,'Color',0.0*[1.0 1.0 1.0])
hold on
plot(s*SE_indep, [2*SE,-2*SE], 'LineWidth', 1.5, 'Color',0.4*[1.0 1.0 1.0] )
xlabel('time lag (minutes)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('autocorrelation', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
xlim([0 s*maxT])
ylim([-0.2 1.0])

print('-depsc2', 'figures/gps_autocorrelation.eps')