clear
drifters = load('sample_data/rho1_drifters_smoothed_interpolated.mat');

%%%% Watching this manual setting!!!
NDrifters = 9;
maxT = [];
for iDrifter = 1:NDrifters
    maxT(end+1) = max(drifters.t{iDrifter});
end
lastTime = min(maxT);

stride = 3;
indices = find(drifters.t{1} >=0 & drifters.t{1} <= lastTime);
reducedIndices = indices(1:3:length(indices));

t = drifters.t{1}(reducedIndices);
x = zeros(length(t),NDrifters);
y = zeros(length(t),NDrifters);

for iDrifter = 1:NDrifters
    indices = find(drifters.t{iDrifter} >=0 & drifters.t{iDrifter} <= lastTime);
    reducedIndices = indices(1:3:length(indices));
    reducedIndices = reducedIndices(1:length(t));
%     drifters.t{iDrifter}(reducedIndices(1))
    x(:,iDrifter) = drifters.x{iDrifter}(reducedIndices);
    y(:,iDrifter) = drifters.y{iDrifter}(reducedIndices);
end

lat0 = drifters.lat0;
lon0 = drifters.lon0;
f0 = drifters.f0;

save('smoothedGriddedRho1Drifters.mat', 't', 'x', 'y', 'lat0', 'lon0', 'f0');