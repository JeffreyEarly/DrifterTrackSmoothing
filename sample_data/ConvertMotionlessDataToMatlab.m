D = importdata('sample_data/motionless_garmin_epix.txt','\t');

lat = D.data(:,1);
lon = D.data(:,2);
lat0 = min(lat) + ( max(lat)-min(lat) )/2;
lon0 = min(lon) + ( max(lon)-min(lon) )/2;
[x,y] = latlon2xy(lat,lon,lat0,lon0);
x = x*1000; y = y*1000;

t = datetime(D.textdata(3:end,2));

save('motionless_garmin_epix.mat', 'x', 'y', 't')