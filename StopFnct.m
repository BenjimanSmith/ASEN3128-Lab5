% This function is used by ODE45 in options to end the drones
% path after a certain tolerance of values is reached
%
%   Author: Benjiman Smith
%   Collaborators: E. Owen, I. Quezada
%   Date: 1/26/2020
%
function [value, isTerm, direction] = StopFnct(t, F)
value = F(12);
isTerm = 1;
direction = [];
end