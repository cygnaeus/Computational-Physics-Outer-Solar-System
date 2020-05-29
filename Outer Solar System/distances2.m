function [r] = distances2(q)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
r = zeros(5,5);
for i = 1:5
    for j = 1:i-1
            a1 = q(3*(i-1) + 1) - q(3*(j-1) + 1);
            a2 = q(3*(i-1) + 2) - q(3*(j-1) + 2);
            a3 = q(3*(i-1) + 3) - q(3*(j-1) + 3);
            r(i,j) = (a1^2 + a2^2 + a3^2)^(-1/2);
            r(j,i) = (a1^2 + a2^2 + a3^2)^(-1/2);
    end
end