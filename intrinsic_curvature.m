function curv = intrinsic_curvature(t,N)
% Return the pointwise intrinsic curvature along the filament, here set to be
% a helix but readily capable of time and space dependence as well as
% constancy. Remove the last line to enable intrinsic curvature in
% simulations.
	s = (0:1/N:1)';
	% Specify the intrinsic curvatures in the d1,d2,d3 directions.
	d1_comp = 0 * s;
	d2_comp = 0 * s;
	d3_comp = 0 * s;
	% Package the output.
	curv = [d1_comp, d2_comp, d3_comp];
	% Remove this line to include the example intrinsic curvature.
	curv = 0 * curv;
end