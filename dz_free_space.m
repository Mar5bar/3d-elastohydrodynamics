function R = dz_free_space(t,Z,EH,N,epsilon,Rs,clamped)
%% DZ_FREE_SPACE(T,Z,EH,N,EPSILON,Rs,CLAMPED) returns the time derivative of Z
%   for the filament given in Z, with current local rotations Rs. Filament
%   aspect ratio is EPSILON, with number of segments N and elastohydrodynamic
%   number EH. T is the current time. Filament is clamped if CLAMPED = TRUE;
	coder.extrinsic('integrated_internal_moments')
	coder.extrinsic('intrinsic_curvature')
    coder.extrinsic('pagemtimes')
    coder.extrinsic('dLab_from_DLocal')

	epsquared = epsilon^2;
	Nsquared = N^2;

    % Generate the directors di.
    [D1, D2, D3] = directors_DLocal(Z);
    d1 = zeros(3,N);
    d2 = zeros(3,N);
    d3 = zeros(3,N);
    [d1, d2, d3] = dLab_from_DLocal(D1, D2, D3, Rs);
    X = spatial_coords_serial(Z(1:3), d3);

    % We will need the local theta and phi values later.
	THETA = Z(4:3+N)';
    PHI = Z(4+N:2*N+3)';

    % Interpolate the directors at the segment endpoints.
    segmentBases = zeros(3,3,N);
    for i = 1 : N
        segmentBases(:,:,i) = [d1(:,i),d2(:,i),d3(:,i)];
    end

    % Approximate via elementwise averages, rescaling singular values to unity.
    meanRots = zeros(3,3,N+1);
    meanRots(:,:,1) = 1.5 * segmentBases(:,:,1) - 0.5 * segmentBases(:,:,2);
    meanRots(:,:,2:end-1) = 0.5 * (segmentBases(:,:,1:end-1) + segmentBases(:,:,2:end));
    meanRots(:,:,end) = 1.5 * segmentBases(:,:,end) - 0.5 * segmentBases(:,:,end-1);
    for i = 1 : N+1
        [U,~,V] = svd(meanRots(:,:,i));
        meanRots(:,:,i) = U*V';
    end
    d1N = reshape(meanRots(:,1,:),[3,N+1]);
    d2N = reshape(meanRots(:,2,:),[3,N+1]);
    d3N = reshape(meanRots(:,3,:),[3,N+1]);

	%------
	% Generate each of the matrices.
	%------

	% Compute inv(A), linking the linear velocity to the forces.
	Ainv = inv_A_rft(d3N,N,epsquared);
	Atilde = zeros(3*N,N);
	for i = 1 : N
		Atilde(3*(i-1)+1:3*i,i) = d3(:,i);
	end
	% We need to add in the rotation of the filament, in terms of the local
	% rate of rotation. Note we have nondimensionalised by 8*pi*mu, hence the
	% factor of 1/2.
	Ainv = blkdiag(Ainv,epsquared*Atilde/2);

	% We will also need B.
	B = zeros(3*N+3,3*(N+1)+3*N);
	% We first build the 3 force balance equations.
	temp = zeros(3*(N+1),1);
	temp(1:3:end) = 2;
	temp(1) = 1; temp(end-2) = 1;
	B(1,1:3*(N+1)) = temp / (2*Nsquared);
	B(2,1:3*(N+1)) = circshift(B(1,1:3*(N+1)),1);
	B(3,1:3*(N+1)) = circshift(B(2,1:3*(N+1)),1);

    % Intermediate quantity required in constructing the remainder of B.
    inter = (d3(:,2:end) + 2*d3(:,1:end-1)) / (6*Nsquared) + (X(:,2:end-1) + X(:,1:end-2))/(2*N);
    
	% Populate the equations for the moment balance in the d1 direction.
    for i = 1 : N
        i_ind = 3*(i-1) + 4; % Steps of 3, starting at 4.
        
        % Contribution to f_i.
        j = i;
        j_ind = 3*(j-1);
        B(i_ind,j_ind+1:j_ind+3) = cross(d1N(:,i),(X(:,j) - X(:,i))/(2*N) + d3(:,j)/(6*Nsquared));
		
        % Contribution to f_j, j = i+1,...,N
        % Efficiently compute cross products of d1 and inter.
        K = [0,-d1N(3,i),d1N(2,i);d1N(3,i),0,-d1N(1,i);-d1N(2,i),d1N(1,i),0];
        d = zeros(N-i,3);
        d(:,1) = inter(1,i:end) - X(1,i)/N;
        d(:,2) = inter(2,i:end) - X(2,i)/N;
        d(:,3) = inter(3,i:end) - X(3,i)/N;
        crosses = K*d';
        B(i_ind,3*i+1:3*N) = crosses(:);

        % Contribution to f_{N+1}.
        j = N+1;
        j_ind = 3*(j-1);
        B(i_ind,j_ind+1:j_ind+3) = cross(d1N(:,i),(X(:,j-1) - X(:,i))/(2*N) + d3(:,j-1)/(3*Nsquared));

        % Add in the tau contributions.
        j_ind = 3*(N+1)+3*(i-1)+1;
        for j = i : N
            B(i_ind,j_ind) = d1N(1,i)/N;
            B(i_ind,j_ind+1) = d1N(2,i)/N;
            B(i_ind,j_ind+2) = d1N(3,i)/N;
            j_ind = j_ind + 3;
        end
    end
    
    % Populate the equations for the moment balance in the d2 direction.
    for i = 1 : N
		i_ind = 3*(i-1) + 5; % Steps of 3, starting at 5.

		% Contribution to f_i.
        j = i;
        j_ind = 3*(j-1);
        B(i_ind,j_ind+1:j_ind+3) = cross(d2N(:,i),(X(:,j) - X(:,i))/(2*N) + d3(:,j)/(6*Nsquared));
        
        % Contribution to f_j, j = i+1,...,N
        % Efficiently compute cross products of d2 and inter.
        K = [0,-d2N(3,i),d2N(2,i);d2N(3,i),0,-d2N(1,i);-d2N(2,i),d2N(1,i),0];
        d = zeros(N-i,3);
        d(:,1) = inter(1,i:end) - X(1,i)/N;
        d(:,2) = inter(2,i:end) - X(2,i)/N;
        d(:,3) = inter(3,i:end) - X(3,i)/N;
        crosses = K*d';
        B(i_ind,3*i+1:3*N) = crosses(:);

        % Contribution to f_{N+1}.
        j = N+1;
        j_ind = 3*(j-1);
        B(i_ind,j_ind+1:j_ind+3) = cross(d2N(:,i),(X(:,j-1) - X(:,i))/(2*N) + d3(:,j-1)/(3*Nsquared));

        % Add in the tau contributions.
        j_ind = 3*(N+1)+3*(i-1)+1;
        for j = i : N
            B(i_ind,j_ind) = d2N(1,i)/N;
            B(i_ind,j_ind+1) = d2N(2,i)/N;
            B(i_ind,j_ind+2) = d2N(3,i)/N;
            j_ind = j_ind + 3;
        end
    end

    % Populate the equations for the moment balance in the d3 direction.
    for i = 1 : N
		i_ind = 3*(i-1) + 6; % Steps of 3, starting at 6.

		% Contribution to f_i.
        j = i;
        j_ind = 3*(j-1);
        B(i_ind,j_ind+1:j_ind+3) = cross(d3N(:,i),(X(:,j) - X(:,i))/(2*N) + d3(:,j)/(6*Nsquared));
        
        % Contribution to f_j, j = i+1,...,N
        % Efficiently compute cross products of d3 and inter.
        K = [0,-d3N(3,i),d3N(2,i);d3N(3,i),0,-d3N(1,i);-d3N(2,i),d3N(1,i),0];
        d = zeros(N-i,3);
        d(:,1) = inter(1,i:end) - X(1,i)/N;
        d(:,2) = inter(2,i:end) - X(2,i)/N;
        d(:,3) = inter(3,i:end) - X(3,i)/N;
        crosses = K*d';
        B(i_ind,3*i+1:3*N) = crosses(:);

        % Contribution to f_{N+1}.
        j = N+1;
        j_ind = 3*(j-1);
        B(i_ind,j_ind+1:j_ind+3) = cross(d3N(:,i),(X(:,j-1) - X(:,i))/(2*N) + d3(:,j-1)/(3*Nsquared));

        % Add in the tau contributions.
        j_ind = 3*(N+1)+3*(i-1)+1;
        for j = i : N
            B(i_ind,j_ind) = d3N(1,i)/N;
            B(i_ind,j_ind+1) = d3N(2,i)/N;
            B(i_ind,j_ind+2) = d3N(3,i)/N;
            j_ind = j_ind + 3;
        end
    end

	%% We build the matrix Q.

	% Begin with the dot(X) rows.
	Q = zeros(3*(N+1),2*N+3);
    % Fill in the dot(X1) dependence.
    Q(1:3:end,1) = 1;
    Q(2:3:end,2) = 1;
    Q(3:3:end,3) = 1;

    %% Compute intermediate products of the RTs and rate of change vectors.
    theta_comp = [cos(THETA).*cos(PHI); cos(THETA).*sin(PHI); -sin(THETA)];
    phi_comp = [-sin(THETA).*sin(PHI); sin(THETA).*cos(PHI); 0*THETA];
    % Reshape the theta_comp and phi_comp so that multiplication by RTs can be
    % done efficiently using pagemtimes.
    theta_comp = reshape(theta_comp,[3,1,N]);
    phi_comp = reshape(phi_comp,[3,1,N]);

    % Multiply the ith vectors by RTi, the transposed rotation matrix for the ith segment.
    theta_comp = pagemtimes(Rs,'transpose',theta_comp,'none'); theta_comp = reshape(theta_comp,[3,N]);
    phi_comp = pagemtimes(Rs,'transpose',phi_comp,'none'); phi_comp = reshape(phi_comp,[3,N]);

    % Rescale by 1/N.
    theta_comp = theta_comp / N;
    phi_comp = phi_comp / N;

    % Assign theta_comp to Q.
    % Loop over columns.
    for j = 1 : N
        Q(3*j + 1 : 3 : end , 3 + j) = theta_comp(1,j);
        Q(3*j + 2 : 3 : end , 3 + j) = theta_comp(2,j);
        Q(3*j + 3 : 3 : end , 3 + j) = theta_comp(3,j);
    end

    % Assign phi_comp to Q.
    % Loop over columns.
    for j = 1 : N
        Q(3*j + 1 : 3 : end , 3 + N + j) = phi_comp(1,j);
        Q(3*j + 2 : 3 : end , 3 + N + j) = phi_comp(2,j);
        Q(3*j + 3 : 3 : end , 3 + N + j) = phi_comp(3,j);
    end

	% Extend to include psi, which requires no transformation.
	Q = blkdiag(Q,eye(N));
	% We also need cos(THETA)*dPHI/dt. This is the block C.
	Q(3*(N+1)+1:3*(N+1)+N,3+N+1:3+2*N) = diag(cos(THETA));

	% Form the RHS.
	R = zeros(3*N+3,1);
	
    % Compute the derivatives of the directors to compute the elastic
    % restoring force. To do this, in order to suppress artefact oscillations
    % in numerics, we will use a crude approximation of the directors at the
    % nodes, and compute first-order finite differences. This sacrifices
    % inherent symmetry in the numerics, recovered as N -> inf, but greatly
    % improves both performance and overall accuracy.
    d1N = [d1; d1(end,:)];
    d2N = [d2; d2(end,:)];
    d3N = [d3; d3(end,:)];
    d1Nds = [d1N(2,:) - d1N(1,:);d1N(2:end,:) - d1N(1:end-1,:)]*N;
    d2Nds = [d2N(2,:) - d2N(1,:);d2N(2:end,:) - d2N(1:end-1,:)]*N;
    d3Nds = [d3N(2,:) - d3N(1,:);d3N(2:end,:) - d3N(1:end-1,:)]*N;

    % We now compute kappa1,kappa2,kappa3 at the s_i, and subtract off the
    % intrinsic curvature. Now, kappa represents the difference between the
    % curvature and the intrinsic curvature, not simply the curvature.
    intrinsic_curv = zeros(N+1,3);
    intrinsic_curv = intrinsic_curvature(t,N);
    kappa1 = dot(d3N,d2Nds)' - intrinsic_curv(:,1);
    kappa2 = dot(d1N,d3Nds)' - intrinsic_curv(:,2);
    kappa3 = dot(d2N,d1Nds)' - intrinsic_curv(:,3);
    
    % We are assuming moment free at the base.
    kappa1(1) = 0;
    kappa2(1) = 0;
    kappa3(1) = 0;

    sigma = 1;
    R(4:3:end) = kappa1(1:end-1);
    R(5:3:end) = kappa2(1:end-1);
    R(6:3:end) = kappa3(1:end-1)/(1+sigma);

	% Add on the integrated contribution of any internally generated moments.
	internal_moments = zeros(3,N);
	internal_moments = integrated_internal_moments(t,N,d1,d2,d3,d1N,d2N,d3N);
	R(4:end) = R(4:end) - internal_moments(:);

	% Form the linear system.
	lin_sys = -EH*B*Ainv*Q;

	% Add on the background flow contribution, noting that the flow needs to
	% be cast in the current basis.
	back = zeros(3*N+3,1);
	back = EH*B*Ainv*background_flow(X,d3,t);
	R = R - back;

    % If clamped, replace overall force and moment free conditions with
    % dot(x)=dot(y)=dot(z)=dot(theta)=dot(phi)=dot(psi) at the base.
    if clamped
        lin_sys(1:6,:) = 0;
        % dx/dt(0)
        lin_sys(1,1) = 1;
        % dy/dt(0)
        lin_sys(2,2) = 1;
        % dz/dt(0)
        lin_sys(3,3) = 1;
        % dtheta/dt(0)
        lin_sys(4,4) = 1;
        % dphi/dt(0)
        lin_sys(5,4+N) = 1;
        % dpsi/dt(0)
        lin_sys(6,4+2*N) = 1;
        % Set the RHS appropriately to 0.
        R(1:6) = 0;
    end

	% Solve the system.
	R = lin_sys \ R;

end