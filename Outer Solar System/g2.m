function [output] = g2(input, H0, L0)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
output = [Energy(input, 1,1,1)- H0 ; Ang_momentum(input)- L0];
end

