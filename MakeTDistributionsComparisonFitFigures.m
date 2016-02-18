% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

% Drifter to highlight in the final plots
choiceDrifter = 1;


shouldDiscardOutliersInACF = 0;
shouldUseSignedACF = 0;

maxlag = 30;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Large error, large tension case
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% How many data points do we have total
drifters_big = load('smoothed_interpolated_rho1_drifters_PowerOptimizedFinalTension');
load('smoothed_interpolated_rho1_drifters_PowerOptimizedFinalTension');
Ndrifters = length(drifters_big.x);

gaussian_pdf_big = @(z) exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi));
position_pdf_big = @(z) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2));

velocity_pdf_big = @(z) exp(-(z.*z)/(2*a*a))/(a*sqrt(2*pi));
velocity_cdf_big = @(z) 0.5*(1 + erf(z/(a*sqrt(2))));    

ACx_big = zeros(maxlag+1,1);
ACy_big = zeros(maxlag+1,1);
a_big = [];
ax_big = [];
ay_big = [];
n_acf_big = 0;
error_big = [];
dist_error_big = [];
for iDrifter = 1:Ndrifters
    x = drifters_big.x{iDrifter};
    y = drifters_big.y{iDrifter};
    t = drifters_big.t{iDrifter};
    N = length(t);
  
    a_big = [a_big; drifters_big.ax{iDrifter};  drifters_big.ay{iDrifter}];
    ax_big = [ax_big; drifters_big.ax{iDrifter}];
    ay_big = [ay_big; drifters_big.ay{iDrifter}];
    
    x_error = drifters_big.x_error{iDrifter};
    y_error = drifters_big.y_error{iDrifter};

    error_big = [error_big; x_error; y_error];
    dist = sqrt( x_error.*x_error + y_error.*y_error );
    dist_error_big = [dist_error_big; dist];
    
    if iDrifter == choiceDrifter
        rejectedPointIndices = find(dist > outlierCut);
        x_error_choice = x_error;
        y_error_choice = y_error;
        dist_choice = dist;
    end
    
    if shouldUseSignedACF == 1
        x_error = sign(x_error);
        y_error = sign(y_error);
    end
    
    if shouldDiscardOutliersInACF == 1
        x_error(abs(x_error) > 6*sigma) = [];
        y_error(abs(y_error) > 6*sigma) = [];
    end
    
    ACx_big = ACx_big + Autocorrelation(x_error,maxlag);
    ACy_big = ACy_big + Autocorrelation(y_error,maxlag);
    n_acf_big = n_acf_big + length(x_error) + length(y_error);
end

ACx_big = ACx_big/Ndrifters;
ACy_big = ACy_big/Ndrifters;
AC_big = (ACx_big + ACy_big)/2;

% log10(std(a_big)/a)

a = std(a_big);
velocity_pdf_big = @(z) exp(-(z.*z)/(2*a*a))/(a*sqrt(2*pi));
velocity_cdf_big = @(z) 0.5*(1 + erf(z/(a*sqrt(2))));    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Small error, small tension case
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

drifters_small = load('smoothed_interpolated_rho1_drifters_NoTension');
load('smoothed_interpolated_rho1_drifters_NoTension');
Ndrifters = length(drifters_small.x);

gaussian_pdf_small = @(z) exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi));
position_pdf_small = @(z) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2));
w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));

exponential_pdf_small = @(z) exp(-(abs(z))/(a))/(2*a);
velocity_cdf_small = @(z) 0.5*(1 + erf(z/(a*sqrt(2))));    

ACx_small = zeros(maxlag+1,1);
ACy_small = zeros(maxlag+1,1);
a_small = [];
ax_small = [];
ay_small = [];
n_acf_small = 0;
error_small = [];
dist_error_small = [];
for iDrifter = 1:Ndrifters
    x = drifters_small.x{iDrifter};
    y = drifters_small.y{iDrifter};
    t = drifters_small.t{iDrifter};
    N = length(t);

    a_small = [a_small; drifters_small.ax{iDrifter};  drifters_small.ay{iDrifter}];
    ax_small = [ax_small; drifters_small.ax{iDrifter}];
    ay_small = [ay_small; drifters_small.ay{iDrifter}];
    
    x_error = drifters_small.x_error_despiked{iDrifter};
    y_error = drifters_small.y_error_despiked{iDrifter};
    
    error_small = [error_small; x_error; y_error];
    dist_error_small = [dist_error_small; sqrt( x_error.*x_error + y_error.*y_error )];
    
    if shouldUseSignedACF == 1
        x_error = sign(x_error);
        y_error = sign(y_error);
    end
        
    if shouldDiscardOutliersInACF == 1
        x_error(abs(x_error) > 6*sigma) = [];
        y_error(abs(y_error) > 6*sigma) = [];
    end
    
    n_acf_small = n_acf_small + length(x_error) + length(y_error);
    
    ACx_small = ACx_small + Autocorrelation(x_error,maxlag);
    ACy_small = ACy_small + Autocorrelation(y_error,maxlag);
end

ACx_small = ACx_small/Ndrifters;
ACy_small = ACy_small/Ndrifters;
AC_small = (ACx_small + ACy_small)/2;

a = std(a_small);
velocity_pdf_small = @(z) exp(-(z.*z)/(2*a*a))/(a*sqrt(2*pi));

% log10(std(a_small)/a)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Position fit figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

drifters = drifters_big;

figure
subplot(2,2,[1 3])
s = 1/1000;
plot(s*drifters_small.x{choiceDrifter},s*drifters_small.y{choiceDrifter},'b'), hold on
plot(s*drifters_big.x{choiceDrifter},s*drifters_big.y{choiceDrifter},'k')
scatter(s*drifters_big.x_raw{choiceDrifter}(rejectedPointIndices),s*drifters_big.y_raw{choiceDrifter}(rejectedPointIndices),(10)^2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
scatter(s*drifters_big.x_raw{choiceDrifter},s*drifters_big.y_raw{choiceDrifter},5)
xlabel('x (km)')
ylabel('y (km)')

subplot(2,2,2)
plot(drifters_small.t{choiceDrifter}/3600,s*drifters_small.x{choiceDrifter},'b'), hold on
plot(drifters_big.t{choiceDrifter}/3600,s*drifters_big.x{choiceDrifter},'k')
scatter(drifters_big.t_raw{choiceDrifter}(rejectedPointIndices)/3600,s*drifters_big.x_raw{choiceDrifter}(rejectedPointIndices),(10)^2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
scatter(drifters.t_raw{choiceDrifter}/3600,s*drifters.x_raw{choiceDrifter},5)
xlabel('t (hours)')
ylabel('x (km)')

subplot(2,2,4)
plot(drifters_small.t{choiceDrifter}/3600,s*drifters_small.y{choiceDrifter},'b'), hold on
plot(drifters_big.t{choiceDrifter}/3600,s*drifters_big.y{choiceDrifter},'k')
scatter(drifters_big.t_raw{choiceDrifter}(rejectedPointIndices)/3600,s*drifters_big.y_raw{choiceDrifter}(rejectedPointIndices),(10)^2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
scatter(drifters.t_raw{choiceDrifter}/3600,s*drifters.y_raw{choiceDrifter},5)
xlabel('t (hours)')
ylabel('y (km)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Position error histogram
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

subplot(2,2,1)
plot_hist_with_pdf( error_big, position_pdf_big, 60, 100 )
% plot_hist_with_pdf( error_big, gaussian_pdf_big, 60, 100 )

subplot(2,2,2)
plot_hist_with_pdf( error_small, position_pdf_small, 60, 100 )
% plot_hist_with_pdf( error_small, gaussian_pdf_small, 20, 100 )

subplot(2,2,3)
plot_hist_with_pdf( a_big, velocity_pdf_big, 10e-5, 100 )

subplot(2,2,4)
plot_hist_with_pdf( a_small, velocity_pdf_small, 10e-5, 100 )
% plot_hist_with_pdf( a_small, exponential_pdf_small, 10e-5, 50 )


figure
subplot(1,2,1)
plot_hist_with_cdf( a_big, velocity_cdf_big, 10e-5, 50 )

subplot(1,2,2)
plot_hist_with_cdf( a_small, velocity_cdf_small, 10e-5, 50 )

return

x = sort(a_small-mean(a_small));
n = length(x);
y_data = (1:n)'/n;
y = velocity_cdf_small(x);
figure
plot(x,y_data), hold on, plot(x,y)
D = max(abs(y-y_data))

lambda = (sqrt(n) + 0.12 + 0.11/sqrt(n))*D;
j=1:25;
p = sum(2*((-1).^(j-1)).*exp(-2*j.*j*lambda*lambda))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Autocorrelation sequence
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(2,1,1)
plot(AC_small), hold on
plot(AC_big)
subplot(2,1,2)
[p1, Q1] = LjungBoxTest(AC_small, n_acf_small);
[p2, Q2] = LjungBoxTest(AC_big, n_acf_big);
plot(p1), hold on
plot(p2)


