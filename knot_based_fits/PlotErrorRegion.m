function PlotErrorRegion( t, x, epsilon)
e_region = [ x + epsilon; flip(x - epsilon,1) ];
t_region = [ t; flip(t,1) ];
stdColor = [0.8 0.8 0.8];
stdEdgeColor = 'none';
fill(t_region,e_region,stdColor,'EdgeColor',stdEdgeColor)
end