function Ainv = inv_A_rft(X,Y,Z,N,epsquared)
%% INV_A_RFT will return the inverse of the matrix A, as given by RFT. This
%	can doubtless be done in a more efficient manner, but is sufficient here.

	% We assume an aspect ratio of epsilon, and divide by 8*pi*mu to match to
	% non-dimensionalisation.
	CN = -0.5 / (log(2/sqrt(epsquared)) - 0.5);
	CT = -0.25 / (log(2/sqrt(epsquared)) - 0.5);

	Adim = 3*(N+1);
	Ainv = zeros(Adim,Adim);
	I3 = eye(3);
	for i = 1 : N+1
		i_ind = 3*(i-1);
		% Approximate the tangent vector.
		switch i
		case 1
			eT = [X(i+1) - X(i); Y(i+1) - Y(i); Z(i+1) - Z(i)];
		case N+1
			eT = [X(i) - X(i-1); Y(i) - Y(i-1); Z(i) - Z(i-1)];
		otherwise
			eT = [X(i+1) - X(i-1); Y(i+1) - Y(i-1); Z(i+1) - Z(i-1)];
		end
		eT = eT / norm(eT);

		% Generate eT*eT^T.
		M = eT * eT';

		% Populate Ainv.
		Ainv(i_ind+1 : i_ind+3, i_ind+1 : i_ind+3) = -CN * (I3 - M) - CT * M;
	end

end