function [t_knot] = NaturalKnotsForSpline( t, K )

if mod(K,2) == 1
   % Odd spline order, so knots go in between points.
   dt = diff(t);
   
   % This gives us N+1 knot points
   t_knot = [t(1); t(1:end-1)+dt/2; t(end)];
   
   % Now remote start and end knots
   for i=1:((K-1)/2)
      t_knot(2) = [];
      t_knot(end-1) = [];
   end
   
else
    t_knot = t;
    
   % Now remote start and end knots
   for i=1:((K-2)/2)
      t_knot(2) = [];
      t_knot(end-1) = [];
   end
    
end

end