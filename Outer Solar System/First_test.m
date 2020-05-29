clear all
close all


%Initual conditions
m = [1, 9.54786e-4, 2.8558e-4, 4.3e-5, 5.177e-5]'; % masses
Jupiter =   [[-3.5,-3.816,-1.550];          1e-3*[5.655,-4.1,-1.9]];
Saturn =    [[9.0755, -3.045,-1.648];       1e-3*[1.68,4.8,1.9]];
Uranus =    [[8.31,-16.29,-7.25];           1e-3*[3.54,1.37, 0.550]];
Neptune =   [[11.47,-25.729,-10.8169];     1e-3*[2.889,1.145, 0.39677]];



init_state(1:15)  = [0,0,0,     Jupiter(1,:),    Saturn(1,:),    Uranus(1,:),    Neptune(1,:)]; %Positions: Sun, Jupiter, Saturn, Uranus, Neptune
v =                 [0,0,0,     Jupiter(2,:),    Saturn(2,:),    Uranus(2,:),    Neptune(2,:)];      %Velocities
init_state(16:30) = v.*kron(m', [1,1,1]);             %Momenta

time_step = 1e1;
G=2.959e-4;
tol = 1e-14;

time = 0:time_step:1e3*356;
k= length(time);
states = zeros(k, 2*3*5);
states(1,:) = init_state;
states1 = zeros(k, 2*3*5);
states1(1,:) = init_state;

for i =2:k
     states(i,:)     = f_implicit_midpoint(states(i-1,:), m, G, time_step, tol);
%     states(i,:)     = f_gauss6(states1(i-1,:), m, G, time_step, tol);
%     states(i,:) = f_explicit_euler(states(i-1,:), m, G, time_step);    
end

% state = states(end,:);
% for i =2:k
%     state     = f_implicit_midpoint(state, m, G, -time_step, tol);
% end
% init_state-state

er = states-states1;
error = sqrt(sum(er.^2,2));
plot(time, error)


%%
close all 
figure
plot(states(:,1), states(:,2), 'b.');
hold on
for i = 1:4
    plot(states(:,3*i + 1 ) -states(:,1), states(:,3*i + 2)-states(:,2), '.');
end
legend("Sun", "Jupiter", "Saturn", "Uranus","Neptune")
title("Location space")
grid on
axis equal

% %%
% 
% figure
% plot3(states(:,1), states(:,2),states(:,3), '*');
% hold on
% for i = 1:4
%     plot3(states(:,3*i+1) -states(:,1), states(:,3*i + 2)-states(:,2),states(:,3*i + 3)-states(:,3), '.');
% end
% legend("Sun", "Jupiter", "Saturn", "Uranus","Neptune")
% title("Location space")
% grid on
% axis equal
% 
% % 
% % %%
% % close all
% % i=1;
% % comet(states(:,3*i+1), states(:,3*i+2))
% % hold on
% % i=2;
% % comet(states(:,3*i+1), states(:,3*i+2))
