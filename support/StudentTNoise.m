function y = StudentTNoise( sigma, nu, n)

if nu < 2
    binEdges = sigma*linspace(-50,50,2001)';
else
    binEdges = sigma*linspace(-100,100,2001)';
end
binWidths = diff(binEdges);

cdf = tcdf(binEdges/sigma, nu);

[~, ~, bin] = histcounts(rand(n,1),cdf);
bin(bin==0 | bin == length(binWidths+1)) = [];
y = binEdges(bin) + rand(length(bin),1).*binWidths(bin);

end
