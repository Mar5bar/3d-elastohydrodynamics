function [D1,D2,D3] = directors_DLocal(Z)
%% directors_DLocal(Z) returns the directors D using the local parameterisation.
	N = (length(Z)-3)/3;
	THETA = Z(4:3+N);
	PHI = Z(4+N:2*N+3);
	PSI = Z(2*N+4:3*N+3);
	D1 = [-sin(PHI).*cos(PSI)-cos(THETA).*cos(PHI).*sin(PSI),cos(PHI).*cos(PSI)-cos(THETA).*sin(PHI).*sin(PSI),sin(THETA).*sin(PSI)]';
	D2 = [sin(PHI).*sin(PSI)-cos(THETA).*cos(PHI).*cos(PSI),-cos(PHI).*sin(PSI)-cos(THETA).*sin(PHI).*cos(PSI),sin(THETA).*cos(PSI)]';
	D3 = [sin(THETA).*cos(PHI),sin(THETA).*sin(PHI),cos(THETA)]';
end