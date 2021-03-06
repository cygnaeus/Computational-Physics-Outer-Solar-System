function [next_state] = func_explicit_euler(state, m, M, G, time_step)
%UNTITLED4 The classical explicit Euler method.
%   Advances the solution with one time step of size 'time_step'.
%   See the description of G_force for detail on the actual system.

next_state = state + time_step*G_force(state, m, M, G);
end