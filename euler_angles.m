function Z = euler_angles(D1,D2,D3)
%% EULER_ANGLES(D1,D2,D3,N) returns the Euler angles parameterising the
%	directors D1,D2,D3.
	THETA = acos(D3(3,:))'; % Value between 0 and pi.
	PHI = atan2(D3(2,:),D3(1,:))'; % Value between -pi and pi.
	PSI = atan2(D1(3,:),D2(3,:))'; % Value between -pi and pi.

	% Note that this Z requires X(1,:)' prepending.
	Z = [THETA;PHI;PSI];
end