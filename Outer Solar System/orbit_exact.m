function [exact_states] = orbit_exact(times,e, N, tol)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

a = 1;
b = a*sqrt(1-e^2);
Ntime = rem(N*times, 2*pi); %N = L_0/a/b

k = length(times);
E = zeros(k, 1);
E(1) = Newton_angle(0, Ntime(1), e, tol);

for i = 2:k
    %Solve next value of E with Newton's method.
    % Use the previous value of E as starting point.
    E(i) = Newton_angle(E(i-1),Ntime(i), e, tol);
    
    if(E(i) > 2*pi)
        E(i) = E(i) - 2*pi;
    end
end

E_dot = N./(1 - e*cos(E));
cosE = cos(E);
sinE  = sin(E);
exact_states = [a*(cosE-e), b*sinE,...
    - a*sinE.*E_dot, b*cosE.*E_dot];

end

