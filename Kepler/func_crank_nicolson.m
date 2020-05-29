function [next_state] = func_crank_nicolson(state, m, M, G, time_step, tol)
%UNTITLED4 The  Crank-Nicolson method.
%   Advances the solution with one time step of size 'time_step'.
%   See the description of G_force for detail on the actual system.

%   For more details on the method, 
%   see https://en.wikipedia.org/wiki/Crank-Nicolson_method

state0 = state;
state1 = state + time_step*G_force((state + state0)/2, m, M, G);
% Simple fixed point method for the implicit solver
while norm(state0 - state1, 2) > tol
    state0 = state1;
    state1 = state + time_step*0.5...
        *(G_force(state , m, M, G)+G_force(state0, m, M, G));
end
next_state = state1;
end

