clear all
close all
%% lateral
Lambda1 = -2; % 1st lambda value, found from lab calculations
Lambda2 = -20;% 2nd lambda value is 10x Lambda 1, in order to make lambda 1 dominant
Ix = 6.8e-5;   % moment of inertia in the x direction [kg m^2]
Iy = 9.2e-5; % Moment of inertia abt y (kgm^2)
g = 9.81
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
for K3 = (.001:.0001:.08); % loop for differing K3 values
 A= [ -K1/Ix, -K2/Ix, -K2*K3/Ix; 1, 0, 0; 0, g, 0]; % matrix derrived in notes from functional block diagram
 Lambda3(:, i) = eig(A); % take the eigenvector of the matrix to get the eigenvalues
 i = i+1; % itterator
end

figure()
  K3 = (0.001:0.0001:.08); % reinitialize K3 vector
plot(real(Lambda3),imag(Lambda3), 'x'); % plot Eigenvalues 3 vs. K3
 xlabel('Real Axis (1/sec)');
 ylabel('Imaginary Axis (1/sec)');
 title('Locus of Eigenvalues for Liearized Lateral Rotation Dynamics');
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
figure()
plot(real(Lambda6),imag(Lambda6), 'x'); % plot Eigenvalues 3 vs. K3
 xlabel('Real Axis (1/sec)');
 ylabel('Imaginary Axis (1/sec)');
 title('Locus of Eigenvalues for Liearized Longitudional Rotation Dynamics');

 index2 = find(Lambda6(3,:) > -0.8, 1, 'last');
k6 = K6(index2+1);
%%





