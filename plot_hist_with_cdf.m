function plot_hist_with_cdf( data, cdf, histwidth, nbins )


% sorted_data = sort(data);


binwidth = histwidth/nbins;
edges = [-histwidth*10;((-histwidth/2+binwidth):binwidth:(histwidth/2-binwidth))';histwidth*10];
binleft = linspace((-histwidth/2),(histwidth/2-binwidth),nbins)';

xi_left = linspace(-histwidth/2,-histwidth/2+binwidth,10)';
xi_mid = linspace(-histwidth/2+binwidth,histwidth/2-binwidth,100)';
xi_right = linspace(histwidth/2-binwidth,histwidth/2,10)';
xi = [xi_left;xi_mid;xi_right];

count = histcounts(data,edges)';
datacdf = cumtrapz(binleft, count/(length(data)*binwidth));
bar(binleft, datacdf, 'histc'); hold on;
plot(xi,cdf(xi),'LineWidth',2,'Color','magenta')
xlim([min(xi) max(xi)])