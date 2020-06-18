function [x,y,z,PSI] = spatial_coords(Z)
%% SPATIAL_COORDS(Z) computes the spatial representation of the filament in
%	Z. This is with respect to the current basis. The twist angle PSI will
%	also be output.
	N = (length(Z)-3)/3;
	PSI = Z(2*N+4:3*N+3);
	tangents = [sin(Z(4:N+3)).*cos(Z(N+4:2*N+3)),sin(Z(4:N+3)).*sin(Z(N+4:2*N+3)),cos(Z(4:N+3))];
	r = cumsum(tangents)/N + [Z(1),Z(2),Z(3)];
	x = [Z(1);r(:,1)];
	y = [Z(2);r(:,2)];
	z = [Z(3);r(:,3)];
end