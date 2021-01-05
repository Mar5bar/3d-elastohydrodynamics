function ub = background_flow(X,d3,t)
%%BACKGROUND_FLOW computes the background flow at X.

% Currently, no flow is prescribed.
	dim = size(X,2);

	u = zeros(1,dim);
	v = zeros(1,dim);
	w = zeros(1,dim);

	% Omega is the curl of the flow field at segment midpoints.
	omega1 = zeros(1,dim-1);
	omega2 = zeros(1,dim-1);
	omega3 = zeros(1,dim-1);

	% Omega should be halved to give the angular velocity.
	Omega = 0.5*[omega1;omega2;omega3];

	% Arrange the elements in the correct order for output.
	ub = [u;v;w]; ub = ub(:);
	ub = [ub;dot(d3,Omega,1)'];
end