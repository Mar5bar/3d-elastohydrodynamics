function [x,y,z,PSI] = spatial_coords_serial(Z)
%% SPATIAL_COORDS(Z) computes the spatial representation of the filament in
%	Z. This is with respect to the current basis. The twist angle PSI will
%	also be output. It does so in a serial manner so as to be able to be
%	converted to a MEX file.
	N = (length(Z)-3)/3;
	PSI = Z(2*N+4:3*N+3);
	x = zeros(N+1,1); x(1) = Z(1);
	y = zeros(N+1,1); y(1) = Z(2);
	z = zeros(N+1,1); z(1) = Z(3);
	for i = 2 : N+1
		x(i) = x(i-1) + sin(Z(2+i))*cos(Z(2+N+i))/N;
		y(i) = y(i-1) + sin(Z(2+i))*sin(Z(2+N+i))/N;
		z(i) = z(i-1) + cos(Z(2+i))/N;
	end
end