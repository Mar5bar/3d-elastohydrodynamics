try close(95)
catch exception
end
%---------------
% Filament setup.
%---------------

global timeouts

% Number of segments.
N = 50;

% Elastohydrodynamic number.
EH = 150000;

% Aspect ratio.
epsilon = 0.01;

% Threshold for THETA values being too close to 0,pi.
delta = pi/20;

%---------------
% Time settings.
%---------------
% T is the end time of the simulation.
T = 10;
ts = linspace(0,T,1001);

%--------------------
% Initial conditions.
%--------------------

% We will denote the vector of minimal coordinates by Z = (x1,y1,z1,theta1,...,thetaN,phi1,...,phiN).
% For an initially-straight vertical configuration with (x1,y1,z1) = 0, we have Z = 0.
Z = zeros(3*N+3,1);
Z(4:3+N) = pi/6; % Set the theta components.
Z(4+N:2*N+3) = linspace(0,2*pi,N); % Set the phi components.
Z(2*N+4:3*N+3) = 0; % Set the psi components.

% Test points for reorientation.
phis = linspace(0,2*pi,100);
thetas = linspace(0,pi,100);
[thetas, phis] = meshgrid(thetas, phis);
thetas = thetas(:); phis = phis(:);
v1s = sin(thetas).*cos(phis);
v2s = sin(thetas).*sin(phis);
v3s = cos(thetas);

tic
% We will loop through, checking that THETA never gets within delta of 0 or pi.
% If it does, we will select a new basis.
T_achieved = 0;
rot = eye(3); % Initial basis.
% Generate the spatial coordinates.
[x,y,z,PSI] = spatial_coords(Z);
X = zeros(N+1,3,length(ts));
X(:,:,1) = [x,y,z];
last_X = X(:,:,1);
last_Z = Z;
D1 = zeros(N,3,length(ts));
D2 = zeros(N,3,length(ts));
D3 = zeros(N,3,length(ts));
current_time_ind = 1;
% This will enable an initial reorientation, usually a good idea.
flag = true;

det_counter = 0;
timeouts = 0;
while (T_achieved < T) % While we have not finished the simulation.
	
	% Get the current directors.
	[d1,d2,d3] = directors(last_Z);

	Z = last_Z;

	% We check to see if THETA is too small/large:
	if (min(pi-Z(4:N+3)) < delta || min(Z(4:N+3)) < delta || flag)
		det_counter = det_counter + 1;
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

		% Rotate the spatial coordinates to the new basis - we are changing basis.
		new_X = (rot*last_X')'; % Note transposes are for MATLAB multiplication, not maths ops.
		% Transform the directors to the new basis.
		new_d1 = (rot*d1')'; new_d2 = (rot*d2')'; new_d3 = (rot*d3')';
		% We compute the angle parameterisation with respect to this new basis.
		Z = [new_X(1,:)';euler_angles(new_d1,new_d2,new_d3)];
		flag = false;
	end

	% Solve the system in this new coordinate system.
	eventFunc = @(t,Z,varargin) odeabort(t,Z,varargin,N,delta); % Aborts solution if near a singularity.
	progressFunc = @(t,y,flag,varargin) odeprog(t,y,flag,varargin,T_achieved); % Displays a progress bar.
	ode_ops = odeset('OutputFcn',progressFunc,'Events',eventFunc,'Stats','off','AbsTol',1e-5,'RelTol',1e-5);

	% Setup the RHS function. Drastic speedup if using user-compiled mex function.

%---Comment/uncomment these lines to use the compiled MEX function. See
%     README.txt for compilation instructions.
	dZ=@(t,z) dz_free_space(t,z,EH,N,epsilon,rot);
	% dZ=@(t,z) dz_free_space_mex(t,z,EH,N,epsilon,rot);
%---

	% Evaluate the solution with ode15s, 
	sol=ode15s(dZ,[0,T-T_achieved],Z,ode_ops);

	% We now compute the solution at the timepoints requested, retaining Z as returned by the solver.
	old_T_achieved = T_achieved;
	T_achieved = sol.x(end) + old_T_achieved;
	% We generate the spatial coordinates of the exact last step taken by the solver.
	% We need to convert these to the original basis, so multiply by inv(rot)=transpose(rot).
	Z = sol.y(:,end);
	[d1,d2,d3] = directors(Z);
	[x,y,z] = spatial_coords(Z);
	last_X = (transpose(rot)*[x,y,z]')'; % Note ' transposes are for MATLAB multiplication, not maths ops.
	original_d1 = (transpose(rot)*d1')'; original_d2 = (transpose(rot)*d2')'; original_d3 = (transpose(rot)*d3')'; 
	last_Z = [last_X(1,:)';euler_angles(original_d1,original_d2,original_d3)];

	% We use sol.y to evaluate the solution at the timepoints in ts.
	% This will prune the ts such that we are only getting timepoints that we need.
	mask = ts <= T_achieved;
	valid_ts = ts(mask);
    valid_ts = valid_ts(current_time_ind:end)-old_T_achieved;
    valid_ts(1) = 0;
    
    % After this setup, perform the solution evaluation.
	Z_at_ts = deval(sol,valid_ts);
	for i = 1 : length(valid_ts)
		[x,y,z,PSI] = spatial_coords(Z_at_ts(:,i));
		X(:,:,current_time_ind+i-1) = (transpose(rot)*[x,y,z]')';
		[d1,d2,d3] = directors(Z_at_ts(:,i));
        Z = [x(1);y(1);z(1);euler_angles(d1,d2,d3)];
		D1(:,:,current_time_ind+i-1) = (transpose(rot)*d1')';
		D2(:,:,current_time_ind+i-1) = (transpose(rot)*d2')';
		D3(:,:,current_time_ind+i-1) = (transpose(rot)*d3')';
	end

	% Update the current time ind.
	current_time_ind = find(ts<=T_achieved,1,'last');
end

toc

disp(['Number of deterministic rotations performed: ',num2str(det_counter),'.'])
disp(['Number of timeouts: ',num2str(timeouts),'.'])
plot_ans(X,ts,10)