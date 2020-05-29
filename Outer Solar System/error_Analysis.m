clear all
close all


%Initual conditions
m = [1, 9.54786e-4, 2.8558e-4, 4.3e-5, 5.177e-5]'; % masses
Jupiter =   [[-3.5,-3.816,-1.550];          1e-3*[5.655,-4.1,-1.9]];
Saturn =    [[9.0755, -3.045,-1.648];       1e-3*[1.68,4.8,1.9]];
Uranus =    [[8.31,-16.29,-7.25];           1e-3*[3.54,1.37, 0.550]];
Neptune =   [[11.47,-25.729,-10.8169];      1e-3*[2.889,1.145, 0.39677]];



init_state(1:15)  = [0,0,0,     Jupiter(1,:),    Saturn(1,:),    Uranus(1,:),    Neptune(1,:)]; %Positions: Sun, Jupiter, Saturn, Uranus, Neptune
v =                 [0,0,0,     Jupiter(2,:),    Saturn(2,:),    Uranus(2,:),    Neptune(2,:)];      %Velocities
init_state(16:30) = v.*kron(m', [1,1,1]);             %Momenta

time_step = 1e1;
G=2.959e-4;
tol = 1e-10;

r = 10; % ratio of saving the states
time = 0:time_step*r:1e4*365;
k= length(time);
states = zeros(k, 2*3*5);
states(1,:) = init_state;
s = init_state;
states1 = zeros(k, 2*3*5);
states1(1,:) = init_state;
s1 = init_state;


for i = 2:k*r
    s = f_implicit_midpoint(s, m, G, time_step, tol);
    s1 = f_gauss6(s1, m, G, time_step, tol);
    if(rem(i, r) == 1)
        l = fix(i/r);
        states(l+1,:) = s;
        states1(l+1,:) = s1;
    end
end

% for i = 2:k*r
%     s = f_implicit_midpoint(s, m, G, time_step, tol);
%     if(rem(i, r) == 1)
%         l = fix(i/r);
%         states(l+1,:) = s;
%         states1(l+1,:) = f_implicit_midpoint(states1(l,:), m, G, r*time_step, tol);
% 
%     end
% end

% er = (states-states1);
% error = sqrt(sum(er.^2,2));
% plot(time, error)

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
%%
figure
plot(states1(:,1), states1(:,2), 'b.');
hold on
for i = 1:4
    plot(states1(:,3*i + 1 ) -states1(:,1), states1(:,3*i + 2)-states1(:,2), '.');
end
legend("Sun", "Jupiter", "Saturn", "Uranus","Neptune")
title("Location space")
grid on
axis equal

%%
clear all
load("Million years", 'time', 'states', 'time_step', 'G', 'tol', 'm')
k = length(time);

figure
E = zeros(1,k);
r = 30;
for i = 1:r:k
    E(i) = Hamiltonian(states(i,:), m, G);
end

plot(time(1:r:k)/365, (E(1:r:k)-E(1))/E(1), 'r.');
grid on
xlabel('Time $t$ (years)', 'Interpreter', 'latex')
ylabel('Relative error of total energy $H$', 'Interpreter', 'latex')
