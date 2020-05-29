function E = Hamiltonian(state, m, G)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
q = state(1:15);  %generalized positions
p = state(16:30); %generalized momenta

r = distances2(q); %r(i, j) = inverse of distance between q_i and q_j (eulidean norm)

m2 = repelem(1./m, 3);
Ek = 1/2 * p.^2 * m2;
Ep = - 1/2*G*m'*r*m  ; %Factor 1/2 cancels the overcount of the summation.
E = Ek + Ep;
end