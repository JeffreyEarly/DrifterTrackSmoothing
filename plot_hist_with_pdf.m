function plot_hist_with_pdf( data, pdf, histwidth, nbins, color )

if nargin < 5
    color = 'b';
end

binwidth = histwidth/nbins;
edges = [-histwidth*10;((-histwidth/2+binwidth):binwidth:(histwidth/2-binwidth))';histwidth*10];
binleft = linspace((-histwidth/2),(histwidth/2-binwidth),nbins)';

xi_left = linspace(-histwidth/2,-histwidth/2+binwidth,10)';
xi_mid = linspace(-histwidth/2+binwidth,histwidth/2-binwidth,100)';
xi_right = linspace(histwidth/2-binwidth,histwidth/2,10)';
xi = [xi_left;xi_mid;xi_right];


edgedensity = integral(pdf,(histwidth/2-binwidth),2*histwidth)/binwidth;
denfunc = edgedensity*ones(size(xi));
denfunc(11:110) = pdf(xi_mid);
count = histcounts(data,edges)';
g = bar(binleft, count/(length(data)*binwidth), 'histc'); hold on;
g.FaceColor = color;
g.FaceAlpha = 0.5;
plot(xi,denfunc,'LineWidth',2,'Color','magenta')
xlim([min(xi) max(xi)])