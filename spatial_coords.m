function X = spatial_coords(X1,d3)
%% SPATIAL_COORDS(X1,d3) computes the spatial representation of the filament
%with tangent d3 from X1.
	N = size(d3,2);
	X = [X1, cumsum(d3,2)/N + X1];
end