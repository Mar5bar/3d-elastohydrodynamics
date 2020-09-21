function integrated_moments_director = integrated_internal_moments(t,N,d1,d2,d3)
% Return moments integrated from s to 1, which can be specified generically.
% Here, as an example, we prescribe sinusoidal moments in the d1 and d2
% directions. This has been disabled by line 28, which should be removed to
% enable actively generated moments. Integrated moments are returned in the
% local director basis.

	s = (0:1/N:1)';
    % Specify the integrated components in the d1,d2,d3 directions in the form
    % of an antiderivative. For example, if the d2 component of the active
    % moment was cos(s-t), specify d2_comp = @(s) sin(s-t).
	d1_comp = @(s) -1 * cos(s-t);
	d2_comp = @(s) 1 * sin(s-t);
	d3_comp = @(s) 0 * s;

    % Integrate the moments over each segment.
    segment_d1_comp = d1_comp(s(2:end)) - d1_comp(s(1:end-1));
    segment_d2_comp = d2_comp(s(2:end)) - d2_comp(s(1:end-1));
    segment_d3_comp = d3_comp(s(2:end)) - d3_comp(s(1:end-1));

    % Construct the moments integrated over each segment, assuming constant
    % directors on each segment.
    segment_moments = segment_d1_comp .* d1 + segment_d2_comp .* d2 + segment_d3_comp .* d3;

    % Integrate the moments from s=s_i to L.
    integrated_moments = cumsum(segment_moments,'reverse');

    % Return them expressed in the local d1,d2,d3 basis at s=s_i.
    integrated_moments_director = [dot(d1,integrated_moments,2),dot(d2,integrated_moments,2),dot(d3,integrated_moments,2)]';

    % Remove this line to allow internal moments.
    integrated_moments_director(:) = 0;
end