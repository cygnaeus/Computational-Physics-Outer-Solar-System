function [L] = Ang_momentum(state)
%UNTITLED5 Return the angular momentum of the state
%   Detailed explanation goes here
L = state(1)*state(4) - state(2)*state(3);

end

