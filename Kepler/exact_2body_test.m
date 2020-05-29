clear all
close all

e = 0.6;

init_state = [1-e ,0,0,sqrt((1+e)/(1-e))]; %State in format: x, y, x_t, y_t
time_step = 0.05;
tol = 1e-8;

m = 1;
M = 1;
G = 1;

time = 0:time_step:200;
k = length(time);
states = orbit_exact(time, e, 1, 1e-6);

figure
plot(states(1,:), states(2,:), 'b-');
grid on
axis equal

E = zeros(k,1);
for i = 1:k
    E(i) = Energy(states(:,1), m, M, G);
end

figure
plot(time, E, 'b-');
grid on
 