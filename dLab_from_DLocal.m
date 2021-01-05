function [d1, d2, d3] = dLab_from_DLocal(D1,D2,D3,Rs)
%% dLab_from_DLocal computes the rotated directors di from Di and the rotation
%   matrices Rs. Rs has shape (3,3,N).
    N = size(Rs,3);
    temp = zeros(3,1,N);
    temp(:,1,:) = D1;
    d1 = reshape(pagemtimes(Rs,'transpose',temp,'none'),[3,N]);
    temp(:,1,:) = D2;
    d2 = reshape(pagemtimes(Rs,'transpose',temp,'none'),[3,N]);
    temp(:,1,:) = D3;
    d3 = reshape(pagemtimes(Rs,'transpose',temp,'none'),[3,N]);
end