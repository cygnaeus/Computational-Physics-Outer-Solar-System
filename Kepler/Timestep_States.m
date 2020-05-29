clear all
close all

n = 30;
tol = 1e-8;

init_state = [1,0,0,1]; %State in format: x, y, x_t, y_t
time_steps = (logspace(-3.1367, -7, n));
Ek_s = zeros(n,1);
Ep_s = zeros(n,1);

m = 1;
M = 1;
G = 1;

parfor j = 1:n
    time_step = time_steps(j);    
    
    time = 0:time_step:5;
    k= length(time);
    states = zeros(k, 4);
    states(1,:) = init_state;
    Ek = zeros(k,1);
    Ep = zeros(k,1);
    [Ek(1), Ep(1)] = Energy(init_state, m, M, G);

    for i =2:k
        states(i,:) = func_implicit_euler(states(i-1,:), m, M, G, time_step, tol)
        %states(i,:) = func_explicit_euler(states(i-1,:), m, M, G, time_step, tol)
        %states(i,:) = func_implicit_midpoint(states(i-1,:), m, M, G, time_step, tol)
        %states(i,:) = func_crank_nicolson(states(i-1,:), m, M, G, time_step, tol)
        [Ek(i), Ep(i)] = Energy(states(i,:), m, M, G);
    end
    Ek_s(j) = Ek(end);
    Ep_s(j) = Ep(end);
    
    disp(j + "/" + n)
end
%%
close all
[Ek1, Ep1] = Energy(init_state, m, M, G);
E = Ek1 + Ep1;


figure
loglog(time_steps, abs((Ek_s + Ep_s)- E), 'bo-');
hold on
plot(time_steps,2e3*time_steps, 'r-', 'Linewidth',1)
grid on
ylim([1e-4, 40])
xlim([8e-8, 1e-3])
title("Implicit Euler")
ylabel("Absolute error of the total energy E_{tot}")
xlabel("Time step h")
legend("Implicit Euler", "O(h)", 'Location', 'southeast')

grid on