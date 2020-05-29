function [change_in_state] = G_force(state, m, M, gamma)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
Const = -gamma*M/(state(1).^2+state(2).^2).^(3/2);

change_in_state = [state(3), state(4), Const/m * state(1), Const/m * state(1)];

change_in_state(3) = Const/m * state(1);
change_in_state(4) = Const/m * state(2);

end

