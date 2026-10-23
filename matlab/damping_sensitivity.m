
clear all,close all,clc;

% parameters
g = 4;
m1 = 1000 * g;
m2 = 800 / g;
k1 = 1.0e6 * g;
k2 = 0.8e6 / g;
c1 = 1.0e6 * g;
c2 = 8.0e6 / g;

% time and force
t = 0:0.001:5;
f = 5000 * sin(2*pi*g*t);

% run original case
a = [0, 1, 0, 0; -(k1+k2)/m1, -(c1+c2)/m1, k2/m1, c2/m1; 0, 0, 0, 1; k2/m2, c2/m2, -k2/m2, -c2/m2];
b = [0; 1/m1; 0; 0];
c = eye(4);
d = zeros(4,1);

sys = ss(a, b, c, d);
y = lsim(sys, f, t);

x1_orig = y(:,1);
x2_orig = y(:,3);
drift_orig = max(abs(x1_orig - x2_orig));

fprintf('original max drift: %.4f m\n', drift_orig);

% case 1: c1 increased 25%
c1_new = c1 * 1.25;
c2_new = c2;

a = [0, 1, 0, 0; -(k1+k2)/m1, -(c1_new+c2_new)/m1, k2/m1, c2_new/m1; 0, 0, 0, 1; k2/m2, c2_new/m2, -k2/m2, -c2_new/m2];
sys = ss(a, b, c, d);
y = lsim(sys, f, t);

x1_c1 = y(:,1);
x2_c1 = y(:,3);
drift_c1 = max(abs(x1_c1 - x2_c1));
change1 = (drift_c1 - drift_orig) / drift_orig * 100;

fprintf('case 1 max drift: %.4f m, change: %.2f%%\n', drift_c1, change1);

% case 2: c2 decreased 25%
c1_new = c1;
c2_new = c2 * 0.75;

a = [0, 1, 0, 0; -(k1+k2)/m1, -(c1_new+c2_new)/m1, k2/m1, c2_new/m1; 0, 0, 0, 1; k2/m2, c2_new/m2, -k2/m2, -c2_new/m2];
sys = ss(a, b, c, d);
y = lsim(sys, f, t);

x1_c2 = y(:,1);
x2_c2 = y(:,3);
drift_c2 = max(abs(x1_c2 - x2_c2));
change2 = (drift_c2 - drift_orig) / drift_orig * 100;

fprintf('case 2 max drift: %.4f m, change: %.2f%%\n', drift_c2, change2);

% case 3: both increased 25%
c1_new = c1 * 1.25;
c2_new = c2 * 1.25;

a = [0, 1, 0, 0; -(k1+k2)/m1, -(c1_new+c2_new)/m1, k2/m1, c2_new/m1; 0, 0, 0, 1; k2/m2, c2_new/m2, -k2/m2, -c2_new/m2];
sys = ss(a, b, c, d);
y = lsim(sys, f, t);

x1_c3 = y(:,1);
x2_c3 = y(:,3);
drift_c3 = max(abs(x1_c3 - x2_c3));
change3 = (drift_c3 - drift_orig) / drift_orig * 100;

fprintf('case 3 max drift: %.4f m, change: %.2f%%\n', drift_c3, change3);

% case 4: both decreased 25%
c1_new = c1 * 0.75;
c2_new = c2 * 0.75;

a = [0, 1, 0, 0; -(k1+k2)/m1, -(c1_new+c2_new)/m1, k2/m1, c2_new/m1; 0, 0, 0, 1; k2/m2, c2_new/m2, -k2/m2, -c2_new/m2];
sys = ss(a, b, c, d);
y = lsim(sys, f, t);

x1_c4 = y(:,1);
x2_c4 = y(:,3);
drift_c4 = max(abs(x1_c4 - x2_c4));
change4 = (drift_c4 - drift_orig) / drift_orig * 100;

fprintf('case 4 max drift: %.4f m, change: %.2f%%\n', drift_c4, change4);


% plot case 1
figure;
subplot(2,2,1);
plot(t, x1_orig*1000, 'b--'); hold on;
plot(t, x2_orig*1000, 'r--');
plot(t, x1_c1*1000, 'b-');
plot(t, x2_c1*1000, 'r-');
xlabel('time (s)');
ylabel('displacement (mm)');
title(sprintf('case 1: c1+25%% (%.2f%% change)', change1));
legend('orig x1', 'orig x2', 'new x1', 'new x2');
grid on;

% plot case 2
subplot(2,2,2);
plot(t, x1_orig*1000, 'b--'); hold on;
plot(t, x2_orig*1000, 'r--');
plot(t, x1_c2*1000, 'b-');
plot(t, x2_c2*1000, 'r-');
xlabel('time (s)');
ylabel('displacement (mm)');
title(sprintf('case 2: c2-25%% (%.2f%% change)', change2));
legend('orig x1', 'orig x2', 'new x1', 'new x2');
grid on;

% plot case 3
subplot(2,2,3);
plot(t, x1_orig*1000, 'b--'); hold on;
plot(t, x2_orig*1000, 'r--');
plot(t, x1_c3*1000, 'b-');
plot(t, x2_c3*1000, 'r-');
xlabel('time (s)');
ylabel('displacement (mm)');
title(sprintf('case 3: both+25%% (%.2f%% change)', change3));
legend('orig x1', 'orig x2', 'new x1', 'new x2');
grid on;

% plot case 4
subplot(2,2,4);
plot(t, x1_orig*1000, 'b--'); hold on;
plot(t, x2_orig*1000, 'r--');
plot(t, x1_c4*1000, 'b-');
plot(t, x2_c4*1000, 'r-');
xlabel('time (s)');
ylabel('displacement (mm)');
title(sprintf('case 4: both-25%% (%.2f%% change)', change4));
legend('orig x1', 'orig x2', 'new x1', 'new x2');
grid on;

% bar chart
figure;
changes = [change1, change2, change3, change4];
bar(changes);
xlabel('case');
ylabel('change in drift (%)');
title('sensitivity analysis results');
xticklabels({'case 1', 'case 2', 'case 3', 'case 4'});
grid on;