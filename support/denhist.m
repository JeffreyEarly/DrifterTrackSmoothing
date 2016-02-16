function dummy = denhist(data, numbins,sel)
% denhist(data, numbins) shows the histogram of the data
% with the number of bins specified by "numbins".
%
% The height of the histogram is equal to 
% (relative freq)/(bin width)
%
  n = length(data);
  binwidth = (max(data)-min(data))/numbins;
  edg= min(data):binwidth:max(data);
  [count, bin] = histc(data, edg);
  
  h = bar(edg, count./(n*binwidth), 'histc');
  color={[0.2 0.4 0.7], [1 0 0.1]};
%   set(h, 'facecolor', color{sel}); % change the color of the bins
%   set(h, 'edgecolor', color{sel});
  
  dummy=count;