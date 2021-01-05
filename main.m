%---------------
% Filament setup.
%---------------

% Use local parameterisations.
new_method_flag = true;
% Reorient velocities of D3 to prolong time away from poles.
reorient_velocities_flag = true;
if new_method_flag
    disp('Using local parameterisations')
    if reorient_velocities_flag
        disp('Reorienting velocities')
    else
        disp('Not reorienting velocities')
    end
else
    disp('Using a global parameterisation')
end


global timeouts

% Number of segments.
N = 50;

% Elastohydrodynamic number.
EH = 150000;

% Aspect ratio.
epsilon = 0.01;

% Threshold for THETA values being too close to 0,pi.
if new_method_flag
	delta = pi/4;
else
	delta = pi/20;
end

%---------------
% Time settings.
%---------------
% T is the end time of the simulation.
T = 10;
ts = linspace(0,T,1001);

% Matrices for saving results.
d1s = zeros(3,N,length(ts));
d2s = zeros(3,N,length(ts));
d3s = zeros(3,N,length(ts));
X = zeros(3,N+1,length(ts));

% Time after which to abort a simulation.
tlim = 60;

%--------------------
% Initial conditions.
%--------------------

% We will denote the vector of minimal coordinates by Z = (x1,y1,z1,theta1,...,thetaN,phi1,...,phiN).
% For an initially-straight vertical configuration with (x1,y1,z1) = 0, we have Z = 0.
Z = zeros(3*N+3,1);
Z(4:3+N) = pi/6; % Set the theta components.
Z(4+N:2*N+3) = linspace(0,2*pi,N); % Set the phi components.
Z(2*N+4:3*N+3) = 0; % Set the psi components.
% These are specified wrt Euler angles in the lab frame!

% Generate the current directors from the parameterisation. We make use of the
% directors_DLocal function only here, as the initial shape is given in lab
% frame Euler angles, yielding d1,d2,d3 rather than D1,D2,D3 as usual.
[d1,d2,d3] = directors_DLocal(Z);

% Is the filament clamped at the base?
clamped = false;

tic
% We will loop through, checking that THETA never gets within delta of 0 or pi.
% If it does, we will select a new basis.
T_achieved = 0;
current_time_ind = 0;

% Create the matrix in which the local rotation for parameterisation will be
% stored.
Rs = zeros(3, 3, N);
I = eye(3);

if ~new_method_flag
    % OLD METHOD------
    % Test points for reorientation.
    phis = linspace(0,2*pi,100);
    thetas = linspace(0,pi,100);
    [thetas, phis] = meshgrid(thetas, phis);
    thetas = thetas(:); phis = phis(:);
    v1s = sin(thetas).*cos(phis);
    v2s = sin(thetas).*sin(phis);
    v3s = cos(thetas);
end

rot_counter = 0;
timeouts = 0;
while (T_achieved < T) % While we have not finished the simulation.
    
	rot_counter = rot_counter + 1;
    if new_method_flag
        % For each segment, compute a rotation matrix R that takes d3 to [1,0,0]^T.
        axs = d3 + [1;0;0]; axs = axs ./ sqrt(sum(axs.^2, 1));
        for i = 1 : N
            v = axs(:,i);
            K = [0,-v(3),v(2);v(3),0,-v(1);-v(2),v(1),0];
            Rs(:,:,i) = I + 2 * K^2;
        end
    else
        %     OLD METHOD FOR COMPARISON
        % Find the test point which is furthest on the sphere from (theta,phi).
        % The locations on the sphere are conveniently given by d3.
        x = [d3(:,1);-d3(:,1)]; y = [d3(:,2);-d3(:,2)]; z = [d3(:,3);-d3(:,3)];
        d = (x-v1s').^2 + (y-v2s').^2 + (z-v3s').^2;
        [~,ind] = max(min(d,[],1));
        theta = thetas(ind);
        phi = phis(ind);
        v = [v1s(ind), v2s(ind), v3s(ind)];
        v = mean([v;0,0,1]); v = v / norm(v);
        % Rotate via Rodrigue's rotation formula.
        K = [0,-v(3),v(2);v(3),0,-v(1);-v(2),v(1),0];
        rot = eye(3) + 2 * K^2;
        for i = 1 : N
            Rs(:,:,i) = rot;
        end
    end

	% Let's rotate the directors to give D1, D2, D3, and compute the resulting
	% parameterisaton.
	[D1,D2,D3] = DLocal_from_dLab(d1,d2,d3,Rs);
	Z = [Z(1:3); euler_angles(D1,D2,D3)];

    % Now that we have the D and parameterisation, compute d/dt(D3) initially
    % to inform a better choice of R.
    if new_method_flag & reorient_velocities_flag
        %---Comment/uncomment these lines to use the compiled MEX function. See
        %     README.txt for compilation instructions.
        dZ_init = dz_free_space(T_achieved,Z,EH,N,epsilon,Rs,clamped);
        % dZ_init = dz_free_space_mex(T_achieved,Z,EH,N,epsilon,Rs,clamped);
        %---
        THETA = Z(4:3+N)'; dTHETA = dZ_init(4:3+N)';
        PHI = Z(4+N:2*N+3)'; dPHI = dZ_init(4+N:2*N+3)';
        dD3dt = [cos(THETA).*cos(PHI); cos(THETA).*sin(PHI); -sin(THETA)] .* dTHETA + ...
                [-sin(THETA).*sin(PHI); sin(THETA).*cos(PHI); 0*THETA] .* dPHI;
        alphas = -atan2(dD3dt(3,:), dD3dt(2,:));
        Rv = @(alpha) [1,0,0;0,cos(alpha),-sin(alpha);0,sin(alpha),cos(alpha)];
        for i = 1 : N
            Rs(:,:,i) = Rv(alphas(i))*Rs(:,:,i);
        end
        % Now reparameterise using the new Rs.
        [D1,D2,D3] = DLocal_from_dLab(d1,d2,d3,Rs);
        Z = [Z(1:3); euler_angles(D1,D2,D3)];
    end

	% Solve the system in Z.
	tstart = tic;
	% Aborts solution if near a singularity or if elapsed time is too great.
	eventFunc = @(t,Z,varargin) odeabort(t,Z,varargin,N,delta,tstart,tlim);
	progressFunc = @(t,y,flag,varargin) odetpbar(t,y,flag); % Displays a progress bar.
	ode_ops = odeset('OutputFcn',progressFunc,'Events',eventFunc,'Stats','off','AbsTol',1e-5,'RelTol',1e-4);

	% Setup the RHS function. Drastic speedup if using user-compiled mex function.

%---Comment/uncomment these lines to use the compiled MEX function. See
%     README.txt for compilation instructions.
	dZ=@(t,z) dz_free_space(t,z,EH,N,epsilon,Rs,clamped);
	% dZ=@(t,z) dz_free_space_mex(t,z,EH,N,epsilon,Rs,clamped);
%---

	% Evaluate the solution with ode15s, 
	sol=ode15s(dZ,[T_achieved,T],Z,ode_ops);

	% We now compute the solution at the timepoints requested, retaining Z as returned by the solver.
	T_achieved = sol.x(end);
	% We generate the spatial coordinates of the exact last step taken by the solver.
	last_Z = sol.y(:,end);
	[last_D1,last_D2,last_D3] = directors_DLocal(last_Z);
	[last_d1,last_d2,last_d3] = dLab_from_DLocal(last_D1,last_D2,last_D3,Rs);
	last_positions = [last_Z(1:3), cumsum(d3,2)/N + last_Z(1:3)];

	% We use sol.y to evaluate the solution at the timepoints in ts.
	% This will prune the ts such that we are only getting timepoints that we need.
	mask = ts <= T_achieved;
	valid_ts = ts(mask);
    valid_ts = valid_ts(current_time_ind+1:end);
    
    % After this setup, perform the solution evaluation.
    if ~isempty(valid_ts)
		Z_at_ts = deval(sol,valid_ts);
		for i = 1 : length(valid_ts)
			[D1,D2,D3] = directors_DLocal(Z_at_ts(:,i));
			[d1,d2,d3] = dLab_from_DLocal(D1,D2,D3,Rs);
			d1s(:,:,current_time_ind+i) = d1;
			d2s(:,:,current_time_ind+i) = d2;
			d3s(:,:,current_time_ind+i) = d3;
			X(:,:,current_time_ind+i) = spatial_coords(Z_at_ts(1:3,i),d3);
		end
	end

	% Prepare for the selection of new Rs and the next ODE solve.
	Z = last_Z;
	d1 = last_d1;
	d2 = last_d2;
	d3 = last_d3;

	% Update the current time ind.
	current_time_ind = find(ts<=T_achieved,1,'last');
end

toc

disp(['Number of times local rotations computed: ',num2str(rot_counter),'.'])
disp(['Number of timeouts: ',num2str(timeouts),'.'])
plot_ans(X,ts,10)