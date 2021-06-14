function H = magnetic_field(t,X)
% Return the prescribed vectorial magnetic field at time t at a point X =
% (x,y,z). X can be a matrix of coordinates, with each column specifying a
% single point.
    H = zeros(size(X));
    H(:,1) = 0*X(:,1);
    H(:,2) = 0*X(:,2);
    H(:,3) = 10*sin(t) + 0*X(:,3);
end