function Z = euler_angles(d1,d2,d3)
%% EULER_ANGLES(d1,d2,d3,N) returns the Euler angles parameterising the
%	directors d1,d2,d3.
	THETA = acos(d3(:,3)); % Value between 0 and pi.
	PHI = atan2(d3(:,2),d3(:,1)); % Value between -pi and pi.
	PSI = atan2(d1(:,3),d2(:,3)); % Value between -pi and pi.

	% Note that this Z requires X(1,:)' prepending.
	Z = [THETA;PHI;PSI];
end