function X = spatial_coords_serial(X1,d3)
%% SPATIAL_COORDS_SERIAL(X1,d3) computes the spatial representation of the
% filament with tangent d3 from X1. It does so in a serial manner so as to be
% able to be converted to a MEX file.
	N = size(d3,2);
	X = zeros(3,N+1); X(:,1) = X1;
    d3 = d3 / N;
    for i = 2 : N+1
        X(:,i) = X(:,i-1) + d3(:,i-1);
    end
end