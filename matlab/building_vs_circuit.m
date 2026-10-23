

clc;
clear;
close all;

% mechanical system (g=4)
m1 = 4000;
m2 = 200;
k1 = 4e6;
k2 = 2e5;
c1 = 4e6;
c2 = 2e6;

% mechanical state space
a_mech = [0, 1, 0, 0;
          -(k1+k2)/m1, -(c1+c2)/m1, k2/m1, c2/m1;
          0, 0, 0, 1;
          k2/m2, c2/m2, -k2/m2, -c2/m2];

b_mech = [0; 1/m1; 0; 0];
c_mech = eye(4);
d_mech = zeros(4, 1);

sys_mech = ss(a_mech, b_mech, c_mech, d_mech);

%electrical system (g=4)
l1 = 4000;       
l2 = 200;       
r1 = 4e6;        
r2 = 2e6;        
c1_elec = 2.5e-7;     
c2_elec = 5e-6;

% electrical state space
a_elec = [0, 1, 0, 0;
          -(1/(l1*c1_elec)+1/(l1*c2_elec)), -(r1+r2)/l1, 1/(l1*c2_elec), r2/l1;
          0, 0, 0, 1;
          1/(l2*c2_elec), r2/l2, -1/(l2*c2_elec), -r2/l2];

b_elec = [0; 1/l1; 0; 0];
c_elec = eye(4);
d_elec = zeros(4, 1);

sys_elec = ss(a_elec, b_elec, c_elec, d_elec);

% simulation
t = linspace(0, 5, 5000);
f = 5000 * sin(8*pi*t);  % force/voltage input
x0 = [0; 0; 0; 0];       % initial conditions

% simulate both systems
y_mech = lsim(sys_mech, f, t, x0);
y_elec = lsim(sys_elec, f, t, x0);

% extract mechanical outputs
x1 = y_mech(:, 1);   % displacement floor 1
v1 = y_mech(:, 2);   % velocity floor 1
x2 = y_mech(:, 3);   % displacement floor 2
v2 = y_mech(:, 4);   % velocity floor 2

% extract electrical outputs
q1 = y_elec(:, 1);   % charge c1 (analogous to displacement)
i1 = y_elec(:, 2);   % current l1 (analogous to velocity)
q2 = y_elec(:, 3);   % charge c2
i2 = y_elec(:, 4);   % current l2
% mechanical velocity = electrical current
err_v1_i1 = abs(v1 - i1);
err_v2_i2 = abs(v2 - i2);

max_err_v1 = max(err_v1_i1);
max_err_v2 = max(err_v2_i2);

percent_v1 = (max_err_v1 / max(abs(v1))) * 100;
percent_v2 = (max_err_v2 / max(abs(v2))) * 100;

fprintf('velocity vs current comparison:\n');
fprintf('v1 vs i1 max error: %.6f (%.3f%%)\n', max_err_v1, percent_v1);
fprintf('v2 vs i2 max error: %.6f (%.3f%%)\n\n', max_err_v2, percent_v2);

% compare displacements with charges
err_x1_q1 = abs(x1 - q1);
err_x2_q2 = abs(x2 - q2);

max_err_x1 = max(err_x1_q1);
max_err_x2 = max(err_x2_q2);

percent_x1 = (max_err_x1 / max(abs(x1))) * 100;
percent_x2 = (max_err_x2 / max(abs(x2))) * 100;

fprintf('displacement vs charge comparison:\n');
fprintf('x1 vs q1 max error: %.6f (%.3f%%)\n', max_err_x1, percent_x1);
fprintf('x2 vs q2 max error: %.6f (%.3f%%)\n\n', max_err_x2, percent_x2);

%% plots - velocity vs current
figure('Position', [100 100 1400 800]);

subplot(2,2,1);
plot(t, v1*1000, 'b-', 'LineWidth', 2); hold on;
plot(t, i1*1000, 'r--', 'LineWidth', 1.5);
xlabel('time (s)');
ylabel('velocity / current (mm/s or mA)');
title(sprintf('floor 1: v1 vs i1 (error: %.3f%%)', percent_v1));
legend('mechanical v1', 'electrical i1');
grid on;

subplot(2,2,2);
plot(t, v2*1000, 'b-', 'LineWidth', 2); hold on;
plot(t, i2*1000, 'r--', 'LineWidth', 1.5);
xlabel('time (s)');
ylabel('velocity / current (mm/s or mA)');
title(sprintf('floor 2: v2 vs i2 (error: %.3f%%)', percent_v2));
legend('mechanical v2', 'electrical i2');
grid on;

subplot(2,2,3);
plot(t, err_v1_i1*1000, 'b', 'LineWidth', 1.5);
xlabel('time (s)');
ylabel('error (mm/s or mA)');
title('v1 vs i1 error');
grid on;

subplot(2,2,4);
plot(t, err_v2_i2*1000, 'r', 'LineWidth', 1.5);
xlabel('time (s)');
ylabel('error (mm/s or mA)');
title('v2 vs i2 error');
grid on;

sgtitle('velocity (mechanical) vs current (electrical) comparison', 'FontWeight', 'bold');

% plots - displacement vs charge
figure('Position', [150 150 1400 800]);

subplot(2,2,1);
plot(t, x1*1000, 'b-', 'LineWidth', 2); hold on;
plot(t, q1*1000, 'r--', 'LineWidth', 1.5);
xlabel('time (s)');
ylabel('displacement / charge (mm or mC)');
title(sprintf('floor 1: x1 vs q1 (error: %.3f%%)', percent_x1));
legend('mechanical x1', 'electrical q1');
grid on;

subplot(2,2,2);
plot(t, x2*1000, 'b-', 'LineWidth', 2); hold on;
plot(t, q2*1000, 'r--', 'LineWidth', 1.5);
xlabel('time (s)');
ylabel('displacement / charge (mm or mC)');
title(sprintf('floor 2: x2 vs q2 (error: %.3f%%)', percent_x2));
legend('mechanical x2', 'electrical q2');
grid on;

subplot(2,2,3);
plot(t, err_x1_q1*1000, 'b', 'LineWidth', 1.5);
xlabel('time (s)');
ylabel('error (mm or mC)');
title('x1 vs q1 error');
grid on;

subplot(2,2,4);
plot(t, err_x2_q2*1000, 'r', 'LineWidth', 1.5);
xlabel('time (s)');
ylabel('error (mm or mC)');
title('x2 vs q2 error');
grid on;

sgtitle('displacement (mechanical) vs charge (electrical) comparison', 'FontWeight', 'bold');

% inter-story drift comparison
drift_mech = abs(x1 - x2);
drift_elec = abs(q1 - q2);

max_drift_mech = max(drift_mech);
max_drift_elec = max(drift_elec);

fprintf('\ninter-story drift:\n');
fprintf('mechanical: %.4f m (%.2f mm)\n', max_drift_mech, max_drift_mech*1000);
fprintf('electrical: %.4f C (%.2f mC)\n', max_drift_elec, max_drift_elec*1000);
fprintf('difference: %.3f%%\n', abs(max_drift_mech-max_drift_elec)/max_drift_mech*100);

figure;
plot(t, drift_mech*1000, 'b-', 'LineWidth', 2); hold on;
plot(t, drift_elec*1000, 'r--', 'LineWidth', 1.5);
xlabel('time (s)');
ylabel('drift (mm or mC)');
title('inter-story drift: mechanical vs electrical');
legend('mechanical |x1-x2|', 'electrical |q1-q2|');
grid on;
