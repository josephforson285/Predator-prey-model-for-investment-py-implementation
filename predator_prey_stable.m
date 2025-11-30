% Parameters
s1 = 0.04;
s2 = 0.05;

K1 = 5;
K2 = 10;

m12 = 0.01;
m21 = 0.02;

v1 = 0.03;
v2 = 0.03;

a1 = 0.1;
a2 = 0.3;

b1 = 0.02;
b2 = 1.5;

c1 = 0.1;
c2 = 0.2;

mu = 0.02;

% Initial values
X1 = 2;
X2 = 8;
Y = 20;

% Time step
dt = 0.1;
iterations = 100000;

% Function defs
r1 = @(X1, X2, Y) 1 / (1 + (a2 * X2 + b2 * Y) / X1);
r2 = @(X1, X2, Y) 1 / (1 + (a1 * X1 + b1 * Y) / X2);

% Derivatives function
function [dX1_dt, dX2_dt, dY_dt] = derivatives(X1, X2, Y, r1, r2, s1, s2, K1, K2, m12, m21, v1, v2, c1, c2, mu)
    dX1_dt = s1 * X1 * (1 - X1 / K1) - m12 * X1 * X2 - v1 * X1 * Y * r1(X1, X2, Y);
    dX2_dt = s2 * X2 * (1 - X2 / K2) - m21 * X1 * X2 - v2 * X2 * Y * r2(X1, X2, Y);
    dY_dt = -mu * Y + c1 * v1 * X1 * Y * r1(X1, X2, Y) + c2 * v2 * X2 * Y * r2(X1, X2, Y);
end

% Runge-Kutta 4th order method step
function [X1_new, X2_new, Y_new] = rk4_step(X1, X2, Y, dt, r1, r2, s1, s2, K1, K2, m12, m21, v1, v2, c1, c2, mu)
    [k1_X1, k1_X2, k1_Y] = derivatives(X1, X2, Y, r1, r2, s1, s2, K1, K2, m12, m21, v1, v2, c1, c2, mu);
    [k2_X1, k2_X2, k2_Y] = derivatives(X1 + 0.5 * dt * k1_X1, X2 + 0.5 * dt * k1_X2, Y + 0.5 * dt * k1_Y, r1, r2, s1, s2, K1, K2, m12, m21, v1, v2, c1, c2, mu);
    [k3_X1, k3_X2, k3_Y] = derivatives(X1 + 0.5 * dt * k2_X1, X2 + 0.5 * dt * k2_X2, Y + 0.5 * dt * k2_Y, r1, r2, s1, s2, K1, K2, m12, m21, v1, v2, c1, c2, mu);
    [k4_X1, k4_X2, k4_Y] = derivatives(X1 + dt * k3_X1, X2 + dt * k3_X2, Y + dt * k3_Y, r1, r2, s1, s2, K1, K2, m12, m21, v1, v2, c1, c2, mu);

    X1_new = X1 + (dt / 6) * (k1_X1 + 2 * k2_X1 + 2 * k3_X1 + k4_X1);
    X2_new = X2 + (dt / 6) * (k1_X2 + 2 * k2_X2 + 2 * k3_X2 + k4_X2);
    Y_new = Y + (dt / 6) * (k1_Y + 2 * k2_Y + 2 * k3_Y + k4_Y);
end

% Preallocate results array
results = zeros(iterations + 1, 3);
results(1, :) = [X1, X2, Y];

% Main loop
for i = 1:iterations
    [X1, X2, Y] = rk4_step(X1, X2, Y, dt, r1, r2, s1, s2, K1, K2, m12, m21, v1, v2, c1, c2, mu);
    results(i + 1, :) = [X1, X2, Y];
end

% Time vector
time = linspace(0, iterations * dt, iterations + 1);

% Plot results
figure;
plot(time, results(:, 1), '-o', 'DisplayName', 'X1');
hold on;
plot(time, results(:, 2), '-o', 'DisplayName', 'X2');
plot(time, results(:, 3), '-o', 'DisplayName', 'Y');
xlabel('Time');
ylabel('Values');
title('Runge-Kutta 4th Order Solution of Differential Equations');
legend;
grid on;

% Print results
fprintf('time | X1       | X2       | Y\n');
for i = 1:iterations + 1
    fprintf('%.1f | %8.5f | %8.5f | %8.5f\n', time(i), results(i, 1), results(i, 2), results(i, 3));
end
