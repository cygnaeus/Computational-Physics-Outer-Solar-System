clear all
close all

n = 30;
tol = 1e-8;

e = 0.6;
init_state = [1-e ,0,0,sqrt((1+e)/(1-e))]; %State in format: x, y, x_t, y_t
time_steps = (logspace(-1, -7, n));
Ek_s = zeros(n,1);
Ep_s = zeros(n,1);

m = 1;
M = 1;
G = 1;

parfor j = 1:n
    time_step = time_steps(j);    
    
    time = 0:time_step:25;
    k= length(time);
    states = zeros(k, 4);
    states(1,:) = init_state;

    for i =2:k
%     state_tilde = func_implicit_midpoint(states(i-1,:), m, M, G, time_step, tol);
    state_tilde = func_explicit_euler(states(i-1,:), m, M, G, time_step);
    %states(i,:) = func_implicit_euler(states(i-1,:), m, M, G, time_step, tol);
    %states(i,:) = func_crank_nicolson(states(i-1,:), m, M, G, time_step, tol);
    
%       states(i,:) = project_H(state_tilde, H0, 1e-6);
%    states(i,:) = project_HL(state_tilde, H0, L0, 1e-10);
    states(i,:) = state_tilde;
    end
    Ek_s(j) = Ek(end);
    Ep_s(j) = Ep(end);
    
    disp(j + "/" + n)
end

timesteps_ExE = time_steps;
globalEr_ExE =  state_norms(exact_sol-states);
smoothGlobal_ExE = smoothdata(globalEr_ExE,'gaussian', floor(length(globalEr_ExE)/125*20));
save('GlobalEr_projection_time.mat', 'timesteps_ExE', 'globalEr_ExE', 'smoothGlobal_ExE', '-append')

%%
clear all
close all
load('GlobalEr_time.mat')

%%%%%Figure with actual and smooth global errrors
% figure
% ax = gca;
% hplot = plot( time_ExE, globalEr_ExE, 'r-' , time_Mid, globalEr_Mid, 'm-',...
%     time_ExE, smoothGlobal_Exe, 'k--', ...
%     time_Mid, smoothGlobal_Mid, 'k-.')
% wid = 2
% % hplot(3).LineWidth = wid;
% % hplot(4).LineWidth = wid;
% labels = ["Explicit Euler $\>$\qquad$h = 0.0001$", "Implicit Midpoint\quad$h = 0.01$", 'Explicit Euler\quad$\>\>\>$ (smoothed)', 'Implicit Midpoint$\>$ (smoothed)'];
%  xlim([0, 125])
%  ylim([0, 3.3])
 
%%%%%Figure with only smooth global errrors
figure
ax = gca;
hplot = plot(time_ExE, smoothGlobal_ExE, 'r-', time_ImE, smoothGlobal_ImE, 'b',...
    time_Mid, smoothGlobal_Mid, 'm-', time_NC, smoothGlobal_NC, 'g-');
 xlim([0, 120])
 ylim([0, 2.5])

hold on
labels = ["Explicit Euler $\>$\qquad$h = 10^{-4}$", "Implicit Euler $\>$\qquad$h = 10^{-4}$", ...
    "Implicit Midpoint\quad$h = 10^{-2}$", "Crank-Nicolson\qquad$h = 10^{-2}$"];
    
grid on
% title("Implicit Euler")
ylabel("Global error of the solution (smoothed)",  'Interpreter','latex')
xlabel("Time $t$" ,'Interpreter','latex')
newOrder = [1 2 4 3];
legend([hplot(newOrder)], [labels(newOrder)], 'Interpreter','latex',   'Location', 'northwest')

grid minor
 set(ax,'XMinorTick','on','YMinorTick','on')
ax.GridAlpha = 0.4;
ax.MinorGridAlpha = 0.1;