function [next_state] = f_implicit_midpoint(state, m,G, time_step, tol)
%UNTITLED4 The  implicit Midpoint method.
%   Advances the solution with one time step of size 'time_step'.
%   See the description of G_force for detail on the actual system.

%   For more details on the method, 
%   see https://en.wikipedia.org/wiki/Midpoint_method


state0 = state;
state1 = state + time_step*G_force5((state + state0)/2, m, G);
% Simple fixed point method for the implicit solver
i = 0;
while norm(state0 - state1, 2) > tol
    state0 = state1;
    state1 = state + time_step*G_force5((state + state0)/2, m, G);
    i = i+1;
    if(i == 1e5)
        error("Fixed point method did not manage to converge in 100 000 steps.\n")
    end
end
next_state = state1;
end


