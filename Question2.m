% This script initializes initial conditions for the trajectory of a drone
% in steady flight conditions (hovering) then calls the ODE45 function to
% numerically integrate the position of the drone with respect to the
% initial conditions, along with plotting the results answering Question
% 7
%
%   Author: Benjiman Smith
%   Collaborators: E. Owen, I. Quezada
%   Date: 1/25/2020
clear all
close all
%% lateral
Lambda1 = -2; % 1st lambda value, found from lab calculations
Lambda2 = -20;% 2nd lambda value is 10x Lambda 1, in order to make lambda 1 dominant
Ix = 6.8e-5;   % moment of inertia in the x direction [kg m^2]
Iy = 9.2e-5; % Moment of inertia abt y (kgm^2)
g = 9.81;
syms K1 K2 % symbollically solving for K1 and K2
eqn1 = Lambda1^2 +Lambda1*(K1/Ix) + (K2/Ix) == 0; % 1st eigenvalue solution
eqn2 = Lambda2^2 +Lambda2*(K1/Ix) + (K2/Ix) == 0; % 2nd eigenvalue solution


[A,B] = equationsToMatrix([eqn1, eqn2], [K1, K2]); % convert both the equations and the unknowns into independent vectors

Solution = linsolve(A,B); % solve the system of equations
Solution = double(Solution); % convert symbolic solution vector to vector of doubles
K1 = Solution(1); % bank control
K2 = Solution(2); % bank control

i = 1;
Lambda3 =zeros(3, 491);
for K3 = (.001:.0001:.08) % loop for differing K3 values
 A= [ -K1/Ix, -K2/Ix, -K2*K3/Ix; 1, 0, 0; 0, g, 0]; % matrix derrived in notes from functional block diagram
 Lambda3(:, i) = eig(A); % take the eigenvector of the matrix to get the eigenvalues
 i = i+1; % itterator
end

%figure(6)
hold on
  K3 = (0.001:0.0001:.08); % reinitialize K3 vector
%plot(K3,real(Lambda3(3,:))); % plot Eigenvalues 3 vs. K3

 index1 = find(Lambda3(3,:) > -0.8, 1, 'last');
k3 = K3(index1+1);
%% longitudional
syms K4 K5 % symbollically solving for K3 and K4
eqn3 = Lambda1^2 +Lambda1*(K4/Iy) + (K5/Iy) == 0; % 1st eigenvalue solution
eqn4 = Lambda2^2 +Lambda2*(K4/Iy) + (K5/Iy) == 0; % 1st eigenvalue solution


[C,D] = equationsToMatrix([eqn3, eqn4], [K4, K5]);

ForceVect2 = linsolve(C,D);
ForceVect2 = double(ForceVect2);
K4 = ForceVect2(1); % elevation control
K5 = ForceVect2(2); % elevation control



i = 1;
for K6 = -(.001:.0001:.08) % loop for differing K3 values
 B= [ -K4/Iy, -K5/Iy, -(K5*K6)/Iy; 1, 0, 0; 0, -g, 0]; % matrix derrived in notes from functional block diagram
 Lambda6(:, i) = eig(B); % take the eigenvector of the matrix to get the eigenvalues
 i = i+1; % itterator
end



  K6 = -(.001:.0001:.08); % reinitialize K3 vector
%plot(K6,real(Lambda6(3,:))); % plot Eigenvalues 3 vs. K3


 index2 = find(Lambda6(3,:) > -0.8, 1, 'last');
k6 = K6(index2+1);
%% feedforeward commands for \Delta v^{E}_{R} 

m = 0.068; % mass of the drone [kg]
r = 0.06;  % body to motor distance [m]
k = 0.0024;  % [Nm/N]
rad = r/sqrt(2);  % [m]
g = 9.81; % gravity [m/s^2]
alpha = 2e-6;   % [N/(m/s)^2]
eta = 1e-3;    % [N/(rad/s)^2]
Ix = 6.8e-5;   % moment of inertia in the x direction [kg m^2]
Iy = 9.2e-5;   % moment of inertia in the y direction [kg m^2]
Iz = 1.35e-4;  % moment of inertia in the z direction [kg m^2]
givens = [alpha eta Ix Iy Iz m r k rad g K1 K2 k3 K4 K5 k6]; % givens vector
F = m*g;
tspan = linspace(0,15); % time vector
Pertubations = zeros(1, 3); % No perturbation
TrimForces = ones(1, 4) * m * g / 4; % forces required by each motor to maintain hover
conditions = zeros(1, 12); % initialize conditions vector (very large);
conditions(12) = -1; % set down direction to 1 to make signs correct for plotting
t1 = 0;
t2 = 0;
X = 0;
X2 = 0;
options = odeset('Events', @StopFnct, 'RelTol', 1e-8); % stop function that ends ODE when a tolerance of 1e-8 is met
[t1, X] = ode45(@(t, F)Specs2LB5LCV(t, F, TrimForces, Pertubations, givens), tspan, conditions, options);
figure()
sgtitle('Change in position with feedforeward commands for \Delta v^{E}_{R}');
subplot(1,3,1);
plot(t1, X(:,10),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('position (m)')
title('Change in x position over time')
subplot(1,3,2);
plot(t1, X(:,11),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('position (m)')
title('Change in y position over time')
subplot(1,3,3);
plot(t1, X(:,12),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('position (m)')
title('Change in z position over time')

figure()
plot3(0, 0, -1, '.', 'MarkerSize',35);
hold on
plot3(X(:, 10), X(:, 11), X(:, 12), 'lineWidth',2)
grid on
zlabel('Down Position, [m]')
xlabel('North Position, [m]')
ylabel('East Position, [m]')
title(' position with feedforeward commands for \Delta v^{E}_{R}');
legend('Origin', 'Change in position over 15 sec');
hold off




figure()
sgtitle('Deviation with feedforeward commands for \Delta v^{E}_{R}');
 subplot(4,2,1);
plot(t1, X(:,7),'linewidth', 2);
 
grid on
xlabel('time (s)')
ylabel('\phi (rad)')
title('change in bank over time')

hold off
%
 subplot(4,2,2);
plot(t1, X(:,8),'linewidth', 2);
 
grid on
xlabel('time (s)')
ylabel('\theta (rad)')
title('change in elevation over time')
hold off
%
 subplot(4,2,3);
plot(t1, X(:,4),'linewidth', 2);
 
grid on
xlabel('time (s)')
ylabel('p (rad/s)')
title('change in roll rate over time')

hold off
%
 subplot(4,2,4);
plot(t1, X(:,5),'linewidth', 2);
 
grid on
xlabel('time (s)')
ylabel('q (rad/s)')
title('change in pitch rate over time')

hold off
%
 subplot(4,2,5);
plot(t1, X(:,6),'linewidth', 2);
 
grid on
xlabel('time (s)')
ylabel('r (rad/s)')
title('change in yaw rate over time')

hold off
%
 subplot(4,2,6);
plot(t1, X(:,1),'linewidth', 2);

grid on
xlabel('time (s)')
ylabel('u (m/s)')
title('change in u over time')

hold off
%
 subplot(4,2,7);
plot(t1, X(:,2),'linewidth', 2);
 
grid on
xlabel('time (s)')
ylabel('v (m/s)')
title('change in v over time')

hold off
%
subplot(4,2,8);
plot(t1, X(:,3),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('w (m/s)')
title('change in w over time')
hold off


%% feedforeward commands for \Delta u^{E}_{R}
m = 0.068; % mass of the drone [kg]
r = 0.06;  % body to motor distance [m]
k = 0.0024;  % [Nm/N]
rad = r/sqrt(2);  % [m]
g = 9.81; % gravity [m/s^2]
alpha = 2e-6;   % [N/(m/s)^2]
eta = 1e-3;    % [N/(rad/s)^2]
Ix = 6.8e-5;   % moment of inertia in the x direction [kg m^2]
Iy = 9.2e-5;   % moment of inertia in the y direction [kg m^2]
Iz = 1.35e-4;  % moment of inertia in the z direction [kg m^2]
givens = [alpha eta Ix Iy Iz m r k rad g K1 K2 k3 K4 K5 k6]; % givens vector
F = m*g;
tspan = linspace(0,15); % time vector
Pertubations = zeros(1, 3); % No perturbation
TrimForces = ones(1, 4) * m * g / 4; % forces required by each motor to maintain hover
conditions = zeros(1, 12); % initialize conditions vector (very large);
conditions(12) = -1; % set down direction to 1 to make signs correct for plotting
t1 = 0;
t2 = 0;
X = 0;
X2 = 0;
options = odeset('Events', @StopFnct, 'RelTol', 1e-8); % stop function that ends ODE when a tolerance of 1e-8 is met
[t1, X] = ode45(@(t, F)Specs2LB5LCU(t, F, TrimForces, Pertubations, givens), tspan, conditions, options);
figure()
sgtitle('Velocity with feedforeward commands for \Delta u^{E}_{R}');
subplot(1,3,1);
plot(t1, X(:,10),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('position (m)')
title('Change in x position over time')
subplot(1,3,2);
plot(t1, X(:,11),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('position (m)')
title('Change in y position over time')
subplot(1,3,3);
plot(t1, X(:,12),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('position (m)')
title('Change in z position over time')

figure()
plot3(0, 0, -1, '.', 'MarkerSize',35);
hold on
plot3(X(:, 10), X(:, 11), X(:, 12), 'lineWidth',2)
grid on
zlabel('Down Position, [m]')
xlabel('North Position, [m]')
ylabel('East Position, [m]')
title(' position with feedforeward commands for \Delta u^{E}_{R}');
legend('Origin', 'Change in position over 15 sec');
hold off
figure()


sgtitle('Deviation with feedforeward commands for \Delta u^{E}_{R}');
 subplot(4,2,1);
plot(t1, X(:,7),'linewidth', 2);
 
grid on
xlabel('time (s)')
ylabel('\phi (rad)')
title('change in bank over time')

hold off
%
 subplot(4,2,2);
plot(t1, X(:,8),'linewidth', 2);
 
grid on
xlabel('time (s)')
ylabel('\theta (rad)')
title('change in elevation over time')
hold off
%
 subplot(4,2,3);
plot(t1, X(:,4),'linewidth', 2);
 
grid on
xlabel('time (s)')
ylabel('p (rad/s)')
title('change in roll rate over time')

hold off
%
 subplot(4,2,4);
plot(t1, X(:,5),'linewidth', 2);
 
grid on
xlabel('time (s)')
ylabel('q (rad/s)')
title('change in pitch rate over time')

hold off
%
 subplot(4,2,5);
plot(t1, X(:,6),'linewidth', 2);
 
grid on
xlabel('time (s)')
ylabel('r (rad/s)')
title('change in yaw rate over time')

hold off
%
 subplot(4,2,6);
plot(t1, X(:,1),'linewidth', 2);

grid on
xlabel('time (s)')
ylabel('u (m/s)')
title('change in u over time')

hold off
%
 subplot(4,2,7);
plot(t1, X(:,2),'linewidth', 2);
 
grid on
xlabel('time (s)')
ylabel('v (m/s)')
title('change in v over time')

hold off
%
subplot(4,2,8);
plot(t1, X(:,3),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('w (m/s)')
title('change in w over time')
hold off


