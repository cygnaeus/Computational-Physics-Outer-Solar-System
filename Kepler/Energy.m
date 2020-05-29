function [E, Ek, Ep] = Energy(state, m, M, G)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Ek = 1/2*m*(state(3).^2 + state(4).^2);
Ep = - G*m*M/sqrt(state(1).^2 + state(2).^2);
E = Ek + Ep;
end

