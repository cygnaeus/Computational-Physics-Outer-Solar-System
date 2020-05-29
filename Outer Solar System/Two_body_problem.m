% clear all
close all

e = 0.6;

init_state = [1-e ,0,0,sqrt((1+e)/(1-e))]; %State in format: x, y, x_t, y_t
time_step = 5e-3;
tol = 1e-14;

m = 1;
M = 1;
G = 1;

time = 0:time_step:(100*2*pi); %200
k= length(time);
states = zeros(k, 4);
states(1,:) = init_state;
exact_sol = orbit_exact(time, e, 1, 1e-14);

Ek = zeros(k,1);
Ep = zeros(k,1);
E  = zeros(k,1);
[E(1), Ek(1), Ep(1)] = Energy(init_state, m, M, G);

for i =2:k
    states(i,:) = func_implicit_midpoint(states(i-1,:), m, M, G, time_step, tol);
%     states(i,:) = func_explicit_euler(states(i-1,:), m, M, G, time_step);
%     states(i,:) = func_implicit_euler(states(i-1,:), m, M, G, time_step, tol);
%     states(i,:) = func_crank_nicolson(states(i-1,:), m, M, G, time_step, tol);
    
    [E(i), Ek(i), Ep(i)] = Energy(states(i,:), m, M, G);
end

% time_NC = time;
% Error_H_NC =  E- E(1);
% state_error_NC =  state_norms(exact_sol-states);
% save('Conservation_of_Energy_time.mat', 'Error_H_NC', 'time_NC', 'state_error_NC', '-append')

% time_ExE = time;
% globalEr_ExE =  state_norms(exact_sol-states);
% smoothGlobal_ExE = smoothdata(globalEr_ExE,'gaussian', floor(length(globalEr_ExE)/125*20));
% save('GlobalEr_time.mat', 'time_ExE', 'globalEr_ExE', 'smoothGlobal_ExE', '-append')

time_Mid = time;
state_Mid = states;
exact_sol2 = exact_sol;
save('Orditals.mat', 'time_Mid', 'state_Mid', 'exact_sol2', '-append')

%%
close all 
% figure
% plot(states(:,1), states(:,2), 'b-');
% title("Location space")
% grid on
% axis equal
% 
% 
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
plot(exact_sol(:,1), exact_sol(:,2), 'r-', 'Linewidth',4);
hold on
plot(states(:,1), states(:,2), 'b-');
grid on
legend("Exact", "Simulation")
subplot(1,2,2)
plot(time, state_norms(exact_sol-states));
hold on 

%plot(time, 3e-4.*time.^2)

plot(time, 16e-7.*time)
%ylim([0,3.5])
legend("Absolute error", "O(t)")
sgtitle("Error evolution - CN method - timestep = " +  time_step);
grid on
% 
% %Max so far of the error
% w = state_norms(exact_sol-states);
% u = zeros(k,1);
% m = zeros(k,1);
% for i = 2:k
%     u(i) = max([u(i-1), w(i)]);
%     i1 = max([1, floor(i-k/125*2*pi)]);
%     m(i) =  mean(w(i1:i));
% end
% figure
% plot(time, u,time, m)
% legend("Maximum error so far", "Mean error so far")

% 
% %%
% clear all
% close all
% load('Conservation_of_Energy_time.mat')
% 
% %%
% close all
% figure
% % o1 = loglog(time_steps_ImMid,1.3e2*time_steps_ImMid, 'k-');
% ax = gca;
% % hold on
% % o2 = plot(time_steps_ImMid,10*time_steps_ImMid.^2, 'k--');
% 
% hplot = plot(time_ExE, Error_H_ExE, 'r-', time_ImE, Error_H_ImE, 'b-', ...
%     time_Mid, Error_H_Mid, 'm-', time_NC, Error_H_NC, 'g-')
% % loglog(time_ImE, Error_Hmax_ImE*2, 'bo-',...
% %     time_ImMid, Error_Hmax_ImMid*2, 'mo-', time_NC, Error_Hmax_NC*2, 'go-' , 'Linewidth',1);
% labels = ["Explicit Euler $\>$\qquad$h = 0.0001$", "Implicit Euler $\>$\qquad$h = 0.0001$", ...
%     "Implicit Midpoint\quad$h = 0.02$", "Crank-Nicolson\qquad$h = 0.02$" ];
% 
% grid on
% % title("Implicit Euler")
% ylabel("Absolute error in total energy $\epsilon_H$",  'Interpreter','latex')
% xlabel("Time $t$" ,'Interpreter','latex')
% newOrder = [ 1 4 3 2];
%  legend([hplot(newOrder)], [labels(newOrder)], 'Interpreter','latex',   'Location', 'northwest')
% 
% grid minor
% axis tight
%  set(ax,'XMinorTick','on','YMinorTick','on')
% ax.GridAlpha = 0.4;
% ax.MinorGridAlpha = 0.1;
% 
% xlim([0, 125])
% ylim([-4e-2, 4e-2])
% 
% %%
% clear all
% close all
% load('GlobalEr_time.mat')
% 
% %%%%%Figure with actual and smooth global errrors
% % figure
% % ax = gca;
% % hplot = plot( time_ExE, globalEr_ExE, 'r-' , time_Mid, globalEr_Mid, 'm-',...
% %     time_ExE, smoothGlobal_Exe, 'k--', ...
% %     time_Mid, smoothGlobal_Mid, 'k-.')
% % wid = 2
% % % hplot(3).LineWidth = wid;
% % % hplot(4).LineWidth = wid;
% % labels = ["Explicit Euler $\>$\qquad$h = 0.0001$", "Implicit Midpoint\quad$h = 0.01$", 'Explicit Euler\quad$\>\>\>$ (smoothed)', 'Implicit Midpoint$\>$ (smoothed)'];
% %  xlim([0, 125])
% %  ylim([0, 3.3])
%  
% %%%%%Figure with only smooth global errrors
% figure
% ax = gca;
% hplot = plot(time_ExE, smoothGlobal_ExE, 'r-', time_ImE, smoothGlobal_ImE, 'b',...
%     time_Mid, smoothGlobal_Mid, 'm-', time_NC, smoothGlobal_NC, 'g-');
%  xlim([0, 120])
%  ylim([0, 2.5])
% 
% hold on
% labels = ["Explicit Euler $\>$\qquad$h = 10^{-4}$", "Implicit Euler $\>$\qquad$h = 10^{-4}$", ...
%     "Implicit Midpoint\quad$h = 10^{-2}$", "Crank-Nicolson\qquad$h = 10^{-2}$"];
%     
% grid on
% % title("Implicit Euler")
% ylabel("Global error of the solution (smoothed)",  'Interpreter','latex')
% xlabel("Time $t$" ,'Interpreter','latex')
% newOrder = [1 2 4 3];
% legend([hplot(newOrder)], [labels(newOrder)], 'Interpreter','latex',   'Location', 'northwest')
% 
% grid minor
%  set(ax,'XMinorTick','on','YMinorTick','on')
% ax.GridAlpha = 0.4;
% ax.MinorGridAlpha = 0.1;
% 
%%
clear all
close all
load('Orditals.mat')

y = 1.5;

ax = figure
plot(exact_sol(:,1), exact_sol(:,2), 'k-', 'Linewidth',2);
hold on
plot(state_ExE_H(:,1), state_ExE_H(:,2), 'b');

grid on
axis equal
ylim([-y,y])
xlim([-2.5,1])
ylabel('y')
xlabel('x')
set(gcf,'position',[0,0,400,300])

% txt = {'Implicit Euler';'h = 0.001'; '50 000 steps'};
% txt2 = {'Implicit Midpoint';'h = 0.05'; '4 000 steps'};
% text([-2.3,0,0], [1.2,1.2,0.9], txt2,'FontSize',10)


txt2 = {'Explicit Euler'; 'with H-projection';'h = 0.005'; '40 000 steps'};
text([-2.3, -2.3,0,0], [1.2, 0.9,1.2,0.9], txt2,'FontSize',10)


