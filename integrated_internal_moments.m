function integrated_moments = integrated_internal_moments(t,N)
% Return moments integrated from s to 1, which can be specified generically.
% Here, as an example, we prescribe sinusoidal moments in the d1 and d2
% directions. This has been disabled by line 13, which should be removed to
% enable actively generated moments.
	s = (0:1/N:1)';
	% Specify the integrated components in the d1,d2,d3 directions.
	d1_comp = 15 * (cos(1 - t) - cos(s - t));
	d2_comp = 15 * (sin(1 - t) - sin(s - t));
	d3_comp = 0 * s;
	% Package the output.
    integrated_moments = [d1_comp(1:end-1), d2_comp(1:end-1), d3_comp(1:end-1)]';
    integrated_moments(:) = 0;
end