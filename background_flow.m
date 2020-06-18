function ub = background_flow(x,y,z,d3,t,row_order,rot)
%%BACKGROUND_FLOW computes the background flow in the current basis, via
%	conversion to the original laboratory frame. 

% Currently, no flow is prescribed.

	% Convert to the laboratory frame.
	X = (transpose(rot)*[x,y,z]')';
	u = zeros(length(x),1);
	v = zeros(length(x),1);;
	w = zeros(length(x),1);;

	% Omega is the curl of the flow field.
	omega1 = zeros(length(d3),1);
	omega2 = zeros(length(d3),1);
	omega3 = zeros(length(d3),1);

	% Convert to the new parameterisation frame.
	Omega = (transpose(rot)*[omega1,omega2,omega3]')';
	ub = (rot*[u, v, w]')';

	% Arrange the elements in the correct order for output.
	ub = ub(row_order);
	ub = [ub; dot(d3,Omega,2)];
end