clear all
close all

e = 0.6;
init_state = [1-e ,0,0,sqrt((1+e)/(1-e))]; %State in format: x, y, x_t, y_t
time_step = 0.001;
tol = 1e-8;

m = 1;
M = 1;
G = 1;

r = 2000; % ratio of saving the states and projection

time = 0:time_step*r:(1e2*2*pi); %200
k= length(time);
states = zeros(k, 4);
states(1,:) = init_state;
state_tilde = init_state;
exact_sol = orbit_exact(time, e, 1, 1e-10);


Ek = zeros(k,1);
Ep = zeros(k,1);
E  = zeros(k,1);
[E(1), Ek(1), Ep(1)] = Energy(init_state, m, M, G);
H0 = E(1);

L0 = Ang_momentum(init_state);

for i = 2:k*r
    state_tilde = func_implicit_midpoint(state_tilde, m, M, G, time_step, tol);
%     state_tilde = func_explicit_euler(states(i-1,:), m, M, G, time_step);
    %states(i,:) = func_implicit_euler(states(i-1,:), m, M, G, time_step, tol);
    %states(i,:) = func_crank_nicolson(states(i-1,:), m, M, G, time_step, tol);
    
%       states(i,:) = project_H(state_tilde, H0, 1e-6);
%    states(i,:) = project_HL(state_tilde, H0, L0, 1e-10);
%     states(i,:) = state_tilde;
%     state_tilde = project_HL(state_tilde, H0, L0, 1e-14);


    if(rem(i, r) == 1)
        l = fix(i/r);
        state_tilde = project_H(state_tilde, H0, 1e-14);
        states(l+1,:) = state_tilde;
    end


    
    %[E(1), Ek(i), Ep(i)] = Energy(states(i,:), m, M, G);
end


time_ExE_H = time;
state_ExE_HL = states;
save('Orditals.mat', 'time_ExE_HL', 'state_ExE_HL', '-append')

%%

% figure
% plot(states(:,1), states(:,2), 'b-');
% title("Location space")
% grid on
% axis equal


% figure
% subplot(1,2,1)
% plot(time, E);
% grid on
% legend("Total energy")
% subplot(1,2,2)
% plot(time, Ek, time, -Ep)
% legend("Kinetic energy", "- Potential energy")
% title("Energy evolution");
% grid on


figure
subplot(1,2,1)
plot(exact_sol(:,1), exact_sol(:,2), 'r.', 'Linewidth',4);
hold on
plot(states(:,1), states(:,2), 'b.');
grid on
legend("Exact", "Simulation")
subplot(1,2,2)
plot(time, state_norms(exact_sol-states), 'b.');
hold on 

%plot(time, 3e-4.*time.^2)

%plot(time, 16e-7.*time)
%ylim([0,3.5])
%legend("Absolute error", "O(t)")
sgtitle("Error evolution - Implicit Midpoint method with HL-projection - timestep = " +  time_step);
grid on

% 
% %Max so far of the error
% w = state_norms(exact_sol-states);
% u = zeros(k,1);
% m = zeros(k,1);
% for i = 2:k
%     u(i) = max([u(i-1), w(i)]);
%     i1 = max([1, floor(i-k*2*pi/200)]);
%     m(i) =  mean(w(i1:i));
% end
% figure
% plot(time, u,time, m)
% legend("Maximum error so far", "Mean error so far")
