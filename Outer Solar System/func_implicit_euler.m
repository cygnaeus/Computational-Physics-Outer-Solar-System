function [next_state] = func_implicit_euler(state, m, M, G, time_step, tol)
%UNTITLED4 The classical explicit Euler method.
%   Advances the solution with one time step of size 'time_step'.
%   See the description of G_force for detail on the actual system.

state0 = state;
state1 = state + time_step*G_force(state0, m, M, G);
% Simple fixed point method for the implicit solver
while norm(state0 - state1, 2) > tol
    state0 = state1;
    state1 = state + time_step*G_force(state0, m, M, G);
end
next_state = state1;
end