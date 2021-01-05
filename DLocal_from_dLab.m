function [D1, D2, D3] = DLocal_from_dLab(d1,d2,d3,Rs)
%% DLocal_from_dLab computes the rotated directors Di from di and the rotation
%   matrices Rs. Rs has shape (3,3,N).
    N = size(Rs,3);
    temp = zeros(3,1,N);
    temp(:,1,:) = d1;
    D1 = reshape(pagemtimes(Rs,temp),[3,N]);
    temp(:,1,:) = d2;
    D2 = reshape(pagemtimes(Rs,temp),[3,N]);
    temp(:,1,:) = d3;
    D3 = reshape(pagemtimes(Rs,temp),[3,N]);
end