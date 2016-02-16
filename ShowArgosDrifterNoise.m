sigma = 1;
nu = 1;

gaussian_pdf = @(z) exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi));
studentt_pdf = @(z) gamma((nu+1)/2)./(sqrt(pi*nu)*sigma*gamma(nu/2)*(1+(z.*z)/(nu*sigma*sigma)).^((nu+1)/2));

n = 1000;
binEdges = linspace(-10,10,201)';
binMids = binEdges(2:end)-(binEdges(2)-binEdges(1))/2;
binWidths = diff(binEdges);
p = studentt_pdf(binMids).*binWidths;
sum(p)

[~, edges, bin] = histcounts(rand(n,1),[0; cumsum(p)]);
bin(bin==0 | bin == length(binWidths+1)) = [];
y = binEdges(bin) + rand(length(bin),1).*binWidths(bin);

figure, histogram(y)

return

% http://www.mathworks.com/matlabcentral/newsreader/view_thread/46607
k = 10;
p = [.1 .2 .3 .4];
binEdges = [0 1 3 6 10];
binWdths = diff(binEdges);
[dum, bin] = histc(rand(1,k),[0 cumsum(p)]);
y = binEdges(bin) + rand(1,k).*binWdths(bin);