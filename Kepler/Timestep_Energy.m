clear all
close all

n = 15;
tol = 1e-8;
e = 0.6;

init_state = [1-e ,0,0,sqrt((1+e)/(1-e))]; %State in format: x, y, x_t, y_t
time_steps = (logspace(-1, -6.9, n)); % -1,-7
Error_H_max = zeros(n,1);

m = 1;
M = 1;
G = 1;

E_initial = Energy(init_state, m, M, G);

parfor j = 1:n
    time_step = time_steps(j);    
    
    time = 0:time_step:25;
    k= length(time);
    state_n = init_state;
    E = zeros(k,1);
    E(1) = Energy(init_state, m, M, G);

    for i =2:k
%         state_n = func_implicit_midpoint(state_n, m, M, G, time_step, tol);
%         state_n = func_crank_nicolson(state_n, m, M, G, time_step, tol);
        state_n = func_explicit_euler(state_n, m, M, G, time_step);
%         state_n = func_implicit_euler(state_n, m, M, G, time_step, tol);

         E(i)  = Energy(state_n, m, M, G);

    end
    Error_H = E - E_initial;
    Error_H_max(j) = max(abs(Error_H));

    disp(j + "/" + n)
end

time_steps_ExE = time_steps;
Error_Hmax_ExE = Error_H_max;
save('Conservation_of_Energy.mat', 'Error_Hmax_ExE', 'time_steps_ExE', '-append')

% %%
% close all
% 
% figure
% loglog(time_steps, Error_H_max, 'bo-');
% hold on
% plot(time_steps,0.85e2*time_steps, 'k--')
% grid on
% % ylim([5e-5, 1])
% % xlim([8e-8, 1e-3])
% % title("Implicit Euler")
% ylabel("Absolute error of the total energy E_{tot}")
% xlabel("Time step h")
% legend("Explicit Euler", "O(h)", 'Location', 'southeast')

grid on

%%
clear all
close all
load('Conservation_of_Energy.mat')

%%
close all
figure
o1 = loglog(time_steps_ImMid,1.3e2*time_steps_ImMid, 'k-');
ax = gca;
hold on
o2 = plot(time_steps_ImMid,10*time_steps_ImMid.^2, 'k--');

hplot = loglog(time_steps_ExE, Error_Hmax_ExE*2, 'ro-', time_steps_ImE, Error_Hmax_ImE*2, 'bo-',...
    time_steps_ImMid, Error_Hmax_ImMid*2, 'mo-', time_steps_NC, Error_Hmax_NC*2, 'go-' , 'Linewidth',1);
labels = {"Explicit Euler", "Implicit Euler", "Implicit Midpoint", "Crank-Nicolson" };

grid on
 ylim([1e-13, 10])
% xlim([8e-8, 1e-3])
% title("Implicit Euler")
ylabel("Relative error in total energy $\mathcal{E}_H$",  'Interpreter','latex')
xlabel("Time step $h$" ,'Interpreter','latex')
newOrder = [2 1 4 3];
legend([hplot(newOrder); o1;o2], [labels(newOrder), "$\mathcal{O}(h)$", "$\mathcal{O}(h^2)$"], 'Interpreter','latex',   'Location', 'southeast')

grid on
ax.GridAlpha = 0.4;
ax.MinorGridAlpha = 0.1;
