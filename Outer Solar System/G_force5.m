function [change_in_state] = G_force5(state, m, G)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
q = state(1:15);  %generalized positions
p = state(16:30); %generalized momenta

[r1, r2, r3] = distances(q); %r_ijk(i, j, k) = k-component of seperation vector q_i -q_j


m2 = repelem(1./m, 3);

change_in_state(1:15) = m2'.* p;
change_in_state(16:3:30) = -G*r1*m.*m;
change_in_state(17:3:30) = -G*r2*m.*m;
change_in_state(18:3:30) = -G*r3*m.*m;

end

