function [next_state] = f_gauss4(state, m, G, time_step, tol)
%UNTITLED4 The  implicit Gauss method of order 4.
%   Advances the solution with one time step of size 'time_step'.
%   Also a 4th order Runge-Kutta.
%   See the description of G_force for detail on the actual system.

%   For more details on the method, 
%   see XXXXXXXXXX

a = 1/4*ones(2) + sqrt(3)/6*[[0 -1];[1 0]];
b = repmat(1/2, 1,2);

k = ones(2,30);
k_new = zeros(2,30);

% Simple fixed point method to solve nonlinear equation
i = 0;
while sum((k - k_new).^2, 'all') > tol^2
    k = k_new;
    k_new(1,:) = G_force5(state + time_step*( a(1,1).*k(1,:)+ a(1,2).*k(2,:)) , m, G);
    k_new(2,:) = G_force5(state + time_step*( a(2,1).*k(1,:)+ a(2,2).*k(2,:)) , m, G);
    i = i+1;
    if(i == 1e5)
        error("Fixed point method did not manage to converge in 100 000 steps.\n")
    end
end
next_state = state + time_step * (b*k)  ;
end