clear all; close all; clc;

%% Given
G = 4;
m1 = 1000 * G;        
m2 = 800 / G;        
k1 = 1.0e6 * G;      
k2 = 0.8e6 / G;       
c1 = 1.0e6 * G;      
c2 = 8.0e6 / G;       


%% State-Space for Mechanical System

A_mech = [0,                    1,              0,              0;
          -(k1+k2)/m1,         -(c1+c2)/m1,     k2/m1,          c2/m1;
          0,                    0,              0,              1;
          k2/m2,                c2/m2,         -k2/m2,         -c2/m2];

B_mech = [0; 1/m1; 0; 0];

C_mech = [1, 0, 0, 0;   
          0, 1, 0, 0;    
          0, 0, 1, 0;    
          0, 0, 0, 1];   

D_mech = zeros(4, 1);

% state-space system
sys_mech = ss(A_mech, B_mech, C_mech, D_mech);

%% time and input force
t = 0:0.001:5;  % 5 seconds
F_input = 5000 * sin(2*pi*G*t); 

%% Simulate Mechanical State-Space Model
[y_mech, t_out, x_mech] = lsim(sys_mech, F_input, t);

% Extract outputs
x1_mech = y_mech(:, 1);  % Floor 1 displacement
v1_mech = y_mech(:, 2);  % Floor 1 velocity
x2_mech = y_mech(:, 3);  % Floor 2 displacement
v2_mech = y_mech(:, 4);  % Floor 2 velocity

%% Calculate Inter-Story Drift
inter_story_drift = abs(x1_mech - x2_mech);
max_drift_original = max(inter_story_drift);
fprintf('Maximum Inter-Story Drift (Original): %.4f m\n\n', max_drift_original);
figure('Name', 'Mechanical System Response');
%%Figure
subplot(3,1,1);
plot(t, F_input/1000, 'k', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Force (kN)');
title('Earthquake Excitation Force');
grid on;

subplot(3,1,2);
plot(t, x1_mech, 'b', 'LineWidth', 1.5); hold on;
plot(t, x2_mech, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Displacement (m)');
title('Floor Displacements vs Time');
legend('Floor 1 (x_1)', 'Floor 2 (x_2)');
grid on;

subplot(3,1,3);
plot(t, v1_mech, 'b', 'LineWidth', 1.5); hold on;
plot(t, v2_mech, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Floor Velocities vs Time');
legend('Floor 1 (v_1)', 'Floor 2 (v_2)');
grid on;

%% Inter-Story Drift
figure('Name', 'Inter-Story Drift');
plot(t, inter_story_drift*1000, 'k', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Inter-Story Drift (mm)');
title(sprintf('Inter-Story Drift |x_1 - x_2| (Max = %.2f mm)', max_drift_original*1000));
grid on;