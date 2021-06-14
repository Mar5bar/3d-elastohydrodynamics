function integrated_moments_director = integrated_internal_moments(t,N,d1,d2,d3,d1N,d2N,d3N,rot,magnetic_moments,X)
% Return moments integrated from s to 1, which can be specified generically.
% Here, we assume that each segment has a magnetic moment, oriented along the
% segment tangent, and there is an external magnetic field. This has been
% disabled by line 23, which should be removed to enable actively generated
% moments. Integrated moments are returned in the local director basis.	

    % Evaluate the magnetic field at the centre of each segment. We rotate the
    % magnetic field to the current frame via rot.
    H = (rot*(magnetic_field(t,0.5*(X(1:end-1,:) + X(2:end,:))))')';

    % The integrated moment on each segment is given by ds *
    % cross(magnetic_moments .* d3, magnetic_field). This is in the lab frame.
    segment_moments = cross(magnetic_moments .* d3,H,2) / N;

    % Integrate the moments from s=s_i to L.
    integrated_moments = cumsum(segment_moments,'reverse');

    % Return them expressed in the local d1,d2,d3 basis at nodes s=s_i.
    integrated_moments_director = [dot(d1N(1:end-1,:),integrated_moments,2),dot(d2N(1:end-1,:),integrated_moments,2),dot(d3N(1:end-1,:),integrated_moments,2)]';

    % Remove this line to allow internal moments.
    integrated_moments_director(:) = 0;
end