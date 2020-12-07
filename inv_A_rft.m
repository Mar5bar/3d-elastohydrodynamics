function Ainv = inv_A_rft(d3N,N,epsquared)
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
		eT = d3N(:,i);

		% Generate eT*eT^T.
		M = eT * eT';

		% Populate Ainv.
		Ainv(i_ind+1 : i_ind+3, i_ind+1 : i_ind+3) = -CN * (I3 - M) - CT * M;
	end

end