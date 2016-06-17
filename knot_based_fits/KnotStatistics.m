addpath('../support')

dt_v = 1;
t=(0:dt_v:100)'; % seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%
% sech^2
%%%%%%%%%%%%%%%%%%%%%%%%%%
u_true = 500;
u_width = 10;
t_0 = 50;
path = @(t) u_true*sech((t-t_0)/u_width).^2;
speed = @(t) -2*(u_true/u_width).*tanh((t-t_0)/u_width).*sech((t-t_0)/u_width).^2;
acceleration = @(t) (u_true/u_width^2)*(4 * (tanh((t-t_0)/u_width)).^2 .* (sech((t-t_0)/u_width)).^2 - 2*(sech((t-t_0)/u_width)).^4);

x_true = path(t);
v_true = speed(t);
a_true = acceleration(t);

sigma = 10; % meters

% Create some Gaussian noise
epsilon = sigma*randn(size(x_true));

% Contaminate the signal
x = x_true + epsilon;

N = length(t);

p = @(z) exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi));
w = @(z)(sigma*sigma);

rng(1)

z_threshold = (2.0:0.2:4.0)';
maxS = 2;
totalRepeats = 25;
stats = zeros(maxS+1,length(z_threshold),totalRepeats);

for iThreshold=1:length(z_threshold)
    for S=0:maxS
        Q_error = [];
        for repeat = 1:totalRepeats
            % Create some Gaussian noise, contaminate the signal
            epsilon = sigma*randn(size(x_true));
            x = x_true + epsilon;

            [t_knot2, ~, constraints] = FindStatisticallySignificantKnotRegions(t,x,sigma,z_threshold(iThreshold),w, S);
            [m_x2,Cm_x2,B2] = bspline_fit_no_tension_constrain(t,x,ones(size(x))*sigma,S,t_knot2,w, constraints);
            tq2 = linspace(t(1),t(end),10*length(x))';
            Bq2 = bspline(tq2,t_knot2,S+1);

            Q_error2 = 0;
            for iS=0:S
                if iS==0
                    thefunction = path;
                elseif iS==1
                    thefunction = speed;
                elseif iS==2
                    thefunction = acceleration;
                end
                x_fit2 = squeeze(Bq2(:,:,iS+1))*m_x2;
                mean_x_error2 = sqrt(mean((thefunction(tq2) - x_fit2).^2));
                Q_error2 = Q_error2 + (mean_x_error2 / sqrt(mean((thefunction(tq2)).^2)) - 1);
            end
            stats(S+1,iThreshold,repeat) = Q_error2;
        end
    end
end

stats_median = median(stats,3);
stats_mean = mean(stats,3);
stats_std = std(stats,0,3);