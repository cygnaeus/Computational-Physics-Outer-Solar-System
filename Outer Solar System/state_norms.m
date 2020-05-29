function [norms] = state_norms(states)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
norms = sqrt(sum(states.^2,2));
end

