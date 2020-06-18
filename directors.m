function [d1,d2,d3] = directors(Z)
%% DIRECTORS(Z) returns the directors with respect to the current basis for
%	the filament in Z.
	N = (length(Z)-3)/3;
	THETA = Z(4:3+N);
	PHI = Z(4+N:2*N+3);
	PSI = Z(2*N+4:3*N+3);
	d1 = [-sin(PHI).*cos(PSI)-cos(THETA).*cos(PHI).*sin(PSI),cos(PHI).*cos(PSI)-cos(THETA).*sin(PHI).*sin(PSI),sin(THETA).*sin(PSI)];
	d2 = [sin(PHI).*sin(PSI)-cos(THETA).*cos(PHI).*cos(PSI),-cos(PHI).*sin(PSI)-cos(THETA).*sin(PHI).*cos(PSI),sin(THETA).*cos(PSI)];
	d3 = [sin(THETA).*cos(PHI),sin(THETA).*sin(PHI),cos(THETA)];
end