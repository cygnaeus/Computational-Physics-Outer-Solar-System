function [next_state] = f_gauss6(state, m, G, time_step, tol)
%UNTITLED4 The  implicit Gauss method of order 6.
%   Advances the solution with one time step of size 'time_step'.
%   Also a 6th order Runge-Kutta.
%   See the description of G_force for detail on the actual system.

%   For more details on the method, 
%   see XXXXXXXXXX

b = [ 5/18 4/9 5/18];
a = repmat(0.5, 3,1)* b + sqrt(15)*[[0, -1/15, -1/30]; [1/24 0 -1/24]; [1/30, 1/15 0]];

k = ones(3,30);
k_new = zeros(3,30);

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