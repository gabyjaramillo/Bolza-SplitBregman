function [value,isterminal,direction] = endpoint(~,y,dirichlet)
% Locate the time when u passes attain boundary value in an increasing direction
% and stop integration.
value = y(1)-dirichlet; % detect u = boundary condition
isterminal = 1;         % stop the integration
direction = 1;          % positive direction
