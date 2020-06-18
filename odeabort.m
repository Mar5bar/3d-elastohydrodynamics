function [value,isterminal,direction]=odeabort(t,Z,varargin,N,delta)
% Adapted from:
% Tim Franklin
% Virginia Tech
% Jesse Norris
% Wake Forrest
% May 2006

%Test to see if we've manually aborted by closing progress bar.
value(1)=double(ishandle(95));
isterminal(1) = 1;
direction(1) = 0;

% We test if our values of THETA are too close to 0 or pi.
value(2)=min(pi-Z(4:N+3)) - delta;
isterminal(2) = 1;
direction(2) = 0;

value(3)=min(Z(4:N+3)) - delta;
isterminal(3) = 1;
direction(3) = 0;