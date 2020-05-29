function [r_ij1, r_ij2, r_ij3] = distances(q)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
r_ij1 = zeros(5,5);
r_ij2 = ones(5,5);
r_ij3 = 2*ones(5,5);

% for i = 1:5
%     for j = 1:i-1
%             a1 = q(3*(i-1) + 1) - q(3*(j-1) + 1);
%             a2 = q(3*(i-1) + 2) - q(3*(j-1) + 2);
%             a3 = q(3*(i-1) + 3) - q(3*(j-1) + 3);
%             
%             r = (a1^2 + a2^2 + a3^2)^(-3/2);
%             r_ij1(i,j) = a1*r;
%             r_ij1(j,i) = -a1*r;
%             r_ij2(i,j) = a2*r;
%             r_ij2(j,i) = -a2*r;
%             r_ij3(i,j) = a3*r;
%             r_ij3(j,i) = -a3*r;
%     end
end

