function R = dz_free_space(t,Z,EH,N,epsilon,rot,clamped)
%% DZ_FREE_SPACE(T,Z,EH,N,EPSILON,ROT,CLAMPED) returns the time derivative of
%   Z for the filament given in Z, with current basis ROT. Filament aspect
%   ratio is EPSILON, with number of segments N and elastohydrodynamic number
%   EH. T is the current time. Filament is clamped if CLAMPED = TRUE;
	coder.extrinsic('integrated_internal_moments')
	coder.extrinsic('intrinsic_curvature')

	epsquared = epsilon^2;
	Nsquared = N^2;

	% Generate the spatial coordinates from the vector Z.
	[x,y,z,~] = spatial_coords_serial(Z);
	% X will store these conveniently.
	X = [x,y,z];

	% We are going to need the directors d1, d2, d3, so we form these from
	% theta, phi. Note the directors here are those on each segment. Pointwise
	% evaluation at nodes will follow.
	THETA = Z(4:3+N);
	PHI = Z(4+N:2*N+3);
	PSI = Z(2*N+4:3*N+3);
	[d1,d2,d3] = directors(Z);

	%------
	% Generate each of the matrices.
	%------

	% Compute inv(A), linking the linear velocity to the forces.
	Ainv = inv_A_rft(x,y,z,N,epsquared);
	Atilde = zeros(3*N,N);
	for i = 1 : N
		Atilde(3*(i-1)+1:3*i,i) = d3(i,:)';
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
    inter = (d3(2:end,:) + 2*d3(1:end-1,:)) / (6*Nsquared) + (X(2:end-1,:) + X(1:end-2,:))/(2*N);
    
	% Populate the equations for the moment balance in the d1 direction.
    for i = 1 : N
        i_ind = 3*(i-1) + 4; % Steps of 3, starting at 4.
        
        % Contribution to f_i.
        j = i;
        j_ind = 3*(j-1);
        B(i_ind,j_ind+1:j_ind+3) = cross(d1(i,:),(X(j,:) - X(i,:))/(2*N) + d3(j,:)/(6*Nsquared));
		
        % Contribution to f_j, j = i+1,...,N
        % Efficiently compute cross products of d1 and inter.
        K = [0,-d1(i,3),d1(i,2);d1(i,3),0,-d1(i,1);-d1(i,2),d1(i,1),0];
        d = zeros(N-i,3);
        d(:,1) = inter(i:end,1) - X(i,1)/N;
        d(:,2) = inter(i:end,2) - X(i,2)/N;
        d(:,3) = inter(i:end,3) - X(i,3)/N;
        crosses = K*d';
        B(i_ind,3*i+1:3*N) = crosses(:);

        % Contribution to f_{N+1}.
        j = N+1;
        j_ind = 3*(j-1);
        B(i_ind,j_ind+1:j_ind+3) = cross(d1(i,:),(X(j-1,:) - X(i,:))/(2*N) + d3(j-1,:)/(3*Nsquared));

        % Add in the tau contributions.
        j_ind = 3*(N+1)+3*(i-1)+1;
        for j = i : N
            B(i_ind,j_ind) = d1(i,1)/N;
            B(i_ind,j_ind+1) = d1(i,2)/N;
            B(i_ind,j_ind+2) = d1(i,3)/N;
            j_ind = j_ind + 3;
        end
    end
    
    % Populate the equations for the moment balance in the d2 direction.
    for i = 1 : N
		i_ind = 3*(i-1) + 5; % Steps of 3, starting at 5.

		% Contribution to f_i.
        j = i;
        j_ind = 3*(j-1);
        B(i_ind,j_ind+1:j_ind+3) = cross(d2(i,:),(X(j,:) - X(i,:))/(2*N) + d3(j,:)/(6*Nsquared));
        
        % Contribution to f_j, j = i+1,...,N
        % Efficiently compute cross products of d2 and inter.
        K = [0,-d2(i,3),d2(i,2);d2(i,3),0,-d2(i,1);-d2(i,2),d2(i,1),0];
        d = zeros(N-i,3);
        d(:,1) = inter(i:end,1) - X(i,1)/N;
        d(:,2) = inter(i:end,2) - X(i,2)/N;
        d(:,3) = inter(i:end,3) - X(i,3)/N;
        crosses = K*d';
        B(i_ind,3*i+1:3*N) = crosses(:);

        % Contribution to f_{N+1}.
        j = N+1;
        j_ind = 3*(j-1);
        B(i_ind,j_ind+1:j_ind+3) = cross(d2(i,:),(X(j-1,:) - X(i,:))/(2*N) + d3(j-1,:)/(3*Nsquared));

        % Add in the tau contributions.
        j_ind = 3*(N+1)+3*(i-1)+1;
        for j = i : N
            B(i_ind,j_ind) = d2(i,1)/N;
            B(i_ind,j_ind+1) = d2(i,2)/N;
            B(i_ind,j_ind+2) = d2(i,3)/N;
            j_ind = j_ind + 3;
        end
    end

    % Populate the equations for the moment balance in the d3 direction.
    for i = 1 : N
		i_ind = 3*(i-1) + 6; % Steps of 3, starting at 6.

		% Contribution to f_i.
        j = i;
        j_ind = 3*(j-1);
        B(i_ind,j_ind+1:j_ind+3) = cross(d3(i,:),(X(j,:) - X(i,:))/(2*N) + d3(j,:)/(6*Nsquared));
        
        % Contribution to f_j, j = i+1,...,N
        % Efficiently compute cross products of d3 and inter.
        K = [0,-d3(i,3),d3(i,2);d3(i,3),0,-d3(i,1);-d3(i,2),d3(i,1),0];
        d = zeros(N-i,3);
        d(:,1) = inter(i:end,1) - X(i,1)/N;
        d(:,2) = inter(i:end,2) - X(i,2)/N;
        d(:,3) = inter(i:end,3) - X(i,3)/N;
        crosses = K*d';
        B(i_ind,3*i+1:3*N) = crosses(:);

        % Contribution to f_{N+1}.
        j = N+1;
        j_ind = 3*(j-1);
        B(i_ind,j_ind+1:j_ind+3) = cross(d3(i,:),(X(j-1,:) - X(i,:))/(2*N) + d3(j-1,:)/(3*Nsquared));

        % Add in the tau contributions.
        j_ind = 3*(N+1)+3*(i-1)+1;
        for j = i : N
            B(i_ind,j_ind) = d3(i,1)/N;
            B(i_ind,j_ind+1) = d3(i,2)/N;
            B(i_ind,j_ind+2) = d3(i,3)/N;
            j_ind = j_ind + 3;
        end
    end

	% We build the matrix Q.
	% Begin with the blocks Q11,...,Q33.
	Q = zeros(3*(N+1),2*N+3);
	% These are blocks Qk1.
	Q(1:N+1,1) = 1;
	Q((N+1)+1:2*(N+1),2) = 1;
	Q(2*(N+1)+1:3*(N+1),3) = 1;

	% We build the submatrices.
	Q12 = zeros(N+1,N); Q13 = zeros(N+1,N); Q22 = zeros(N+1,N); Q23 = zeros(N+1,N); Q32 = zeros(N+1,N); Q33 = zeros(N+1,N);
	for j = 1 : N
		Q12(j+1:end,j) =   cos(THETA(j))*cos(PHI(j));
		Q13(j+1:end,j) = - sin(THETA(j))*sin(PHI(j));
		Q22(j+1:end,j) =   cos(THETA(j))*sin(PHI(j));
		Q23(j+1:end,j) =   sin(THETA(j))*cos(PHI(j));
		Q32(j+1:end,j) = - sin(THETA(j));
		% Q33 = 0 identically.
	end
	Q(:,4:end) = [Q12, Q13; Q22, Q23; Q32, Q33] / N;

	row_order = zeros(3*(N+1),1);
	for j = 1 : N + 1
		row_order(3*(j-1)+1:3*j) = [j,N+1+j,2*N+2+j];
	end
	Q = Q(row_order,:);
	% Extend to include psi, which requires no transformation.
	Q = blkdiag(Q,eye(N));
	% We also need cos(THETA)*dPHI/dt. This is the block C.
	Q(3*(N+1)+1:3*(N+1)+N,3+N+1:3+2*N) = diag(cos(THETA));

	% Form the RHS.
	R = zeros(3*N+3,1);
	% 4th order central differences, 2nd order differences at ends.
	d1ds = [-0.5*d1(3,:) + 2*d1(2,:) - 1.5*d1(1,:); 0.5*d1(3,:) - 0.5*d1(1,:); -1/12*d1(5:end,:) + 2/3*d1(4:end-1,:) - 2/3*d1(2:end-3,:) + 1/12*d1(1:end-4,:); 0.5*d1(end,:) - 0.5*d1(end-2,:); 1.5*d1(end,:) - 2*d1(end-1,:) + 0.5*d1(end-2,:)]*N;
	d2ds = [-0.5*d2(3,:) + 2*d2(2,:) - 1.5*d2(1,:); 0.5*d2(3,:) - 0.5*d2(1,:); -1/12*d2(5:end,:) + 2/3*d2(4:end-1,:) - 2/3*d2(2:end-3,:) + 1/12*d2(1:end-4,:); 0.5*d2(end,:) - 0.5*d2(end-2,:); 1.5*d2(end,:) - 2*d2(end-1,:) + 0.5*d2(end-2,:)]*N;
	d3ds = [-0.5*d3(3,:) + 2*d3(2,:) - 1.5*d3(1,:); 0.5*d3(3,:) - 0.5*d3(1,:); -1/12*d3(5:end,:) + 2/3*d3(4:end-1,:) - 2/3*d3(2:end-3,:) + 1/12*d3(1:end-4,:); 0.5*d3(end,:) - 0.5*d3(end-2,:); 1.5*d3(end,:) - 2*d3(end-1,:) + 0.5*d3(end-2,:)]*N;

	% We now compute kappa1,kappa2,kappa3 at the s_i, and subtract off the
	% intrinsic curvature. Now, kappa represents the difference between the
	% curvature and the intrinsic curvature, not simply the curvature.
	intrinsic_curv = zeros(N,3);
	intrinsic_curv = intrinsic_curvature(t,N);
	kappa1 = dot(d3,d2ds,2) - intrinsic_curv(:,1);
	kappa2 = dot(d1,d3ds,2) - intrinsic_curv(:,2);
	kappa3 = dot(d2,d1ds,2) - intrinsic_curv(:,3);
	
	% We are assuming moment free at the base.
	kappa1(1) = 0;
	kappa2(1) = 0;
	kappa3(1) = 0;

	sigma = 1;
	R(4:3:end) = kappa1;
	R(5:3:end) = kappa2;
	R(6:3:end) = kappa3/(1+sigma);

	% Add on the integrated contribution of any internally generated moments.
	internal_moments = zeros(3,N);
	internal_moments = integrated_internal_moments(t,N,d1,d2,d3);
	R(4:end) = R(4:end) - internal_moments(:);

	% Form the linear system.
	lin_sys = -EH*B*Ainv*Q;

	% Add on the background flow contribution, noting that the flow needs to
	% be cast in the current basis.
	back = zeros(3*N+3,1);
	back = EH*B*Ainv*background_flow(x,y,z,d3,t,row_order,rot);
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