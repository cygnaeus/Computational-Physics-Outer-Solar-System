function [y_res] = project_H(y_t, H0, tol);
%UNTITLED2 Find the closest  state with total energy 'H0'
%   Projects the state 'y_t' to the closest state with
%   total energy 'H0'. The distance is defined by the 
%   Euclidean norm.
%
%   The method is to solve the Lagrange multiplier lamda
%   with simplified Newton's iteration. 
%
%   y_t     1x4 state vector to be projected
%   H0      Desired total energy of the output state
%   tol     Absolute tolerance of the Simplified Newton iteration
%   y_res   1x4 output state vector with approximate total energy H0.

x = y_t(1);
y = y_t(2);
x_dot = y_t(3);
y_dot = y_t(4);
r = sqrt(x^2+ y^2);

g1_dot = [x/r^3, y/r^3, x_dot, y_dot];
inv_m =  1/(1/r^4 + x_dot^2+ y_dot^2 );

lamda0 = 0;
lamda1 = lamda0 - inv_m* g1(y_t + g1_dot*lamda0, H0);

while (abs(lamda0 - lamda1) > tol)
    lamda0 = lamda1;
    lamda1 = lamda0 - inv_m* g1(y_t + g1_dot*lamda0, H0);
end

y_res = y_t + g1_dot*lamda1;
end

