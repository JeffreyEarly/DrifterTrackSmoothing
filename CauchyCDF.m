alpha = (0.1:0.1:10)';

alpha = gamma;

cdf_minus_1 = atan( (-1-alpha)./alpha )/pi + 1/2;
cdf_plus_1 = atan( (1-alpha)./alpha )/pi + 1/2;

prob = 1 - (cdf_plus_1 - cdf_minus_1);