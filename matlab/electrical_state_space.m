clc;
clear;
close all;

%G=4
L1 = 4000;        
L2 = 200;        
R1 = 4e6;         
R2 = 2e6;         
C1 = 2.5e-7;      
C2 = 5e-6;        

% State-space matrices

A = [ ...
    0  1  0  0;
   -(1/(L1*C1)+1/(L1*C2))  -(R1+R2)/L1   1/(L1*C2)   R2/L1;
    0  0  0  1;
    1/(L2*C2)   R2/L2   -1/(L2*C2)  -R2/L2 ];

B = [0; 1/L1; 0; 0];
C = eye(4);
D = zeros(4,1);

%state-space system
sys_elec = ss(A,B,C,D);
t = linspace(0,5,5000);              
V = 5000 * sin(8*pi*t);              
% initial conditions
x0 = [0 0 0 0];   % [q1 i1 q2 i2]

% Simulate system

[y, t, x] = lsim(sys_elec, V, t, x0);

% Extract outputs
q1 = y(:,1);      % Charge C1
i1 = y(:,2);      % Current L1
q2 = y(:,3);      % Charge C2
i2 = y(:,4);      % Current L2

% Plot results
figure;
subplot(2,1,1)
plot(t, q1, 'LineWidth', 1.5); hold on;
plot(t, q2, 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Charge (C)');
title('Capacitor Charges (Electrical Displacement)');
legend('q1 (Floor 1)', 'q2 (Floor 2)', 'Interpreter', 'none');

subplot(2,1,2)
plot(t, i1, 'LineWidth', 1.5); hold on;
plot(t, i2, 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Current (A)');
title('Inductor Currents (Electrical Velocity)');
legend('i1 (Floor 1)', 'i2 (Floor 2)', 'Interpreter', 'none');


%Inter-story drift
figure;
plot(t, abs(q1 - q2), 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('|q1 - q2| (C)');
title('Inter-story Drift (Electrical Analog)');
