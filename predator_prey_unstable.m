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
b2 = 3;   
c1 = 0.1;
c2 = 0.23;  
mu = 0.02;

% Initial values
X1_0 = 0.54; 
X2_0 = 3.3; 
Y_0 = 0.75; 

% Initial conditions 
y0 = [X1_0; X2_0; Y_0];

% Time span
tspan = linspace(0, 20000, 20000); 

% function for system of DE
function dydt = prey_predator_system(~, y)
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
    b2 = 3;
    c1 = 0.1;
    c2 = 0.23;
    mu = 0.02;

    X1 = y(1);
    X2 = y(2);
    Y = y(3);

    r1 = X1 / (X1 + a2 * X2 + b2 * Y);
    r2 = X2 / (X2 + a1 * X1 + b1 * Y);

    dX1_dt = s1 * X1 * (1 - X1 / K1) - m12 * X1 * X2 - v1 * X1 * Y * r1;
    dX2_dt = s2 * X2 * (1 - X2 / K2) - m21 * X1 * X2 - v2 * X2 * Y * r2;
    dY_dt = -mu * Y + c1 * v1 * X1 * Y * r1 + c2 * v2 * X2 * Y * r2;

    dydt = [dX1_dt; dX2_dt; dY_dt];
end

% Solve 
[t, y] = ode45(@prey_predator_system, tspan, y0);

% Plot 
figure;
plot(t, y(:,1), 'b-', 'DisplayName', 'X1'); 
hold on;
plot(t, y(:,1), 'bo');
plot(t, y(:,2), 'r-', 'DisplayName', 'X2'); 
plot(t, y(:,2), 'ro'); 
plot(t, y(:,3), 'g-', 'DisplayName', 'Y');
plot(t, y(:,3), 'go'); 
hold off;
xlabel('Time');
ylabel('Values');
title('ODE45 Solution of Differential Equations');
legend;
grid on;

% Save data
data = [t, y];
fileID = fopen('results_1.txt', 'w');
fprintf(fileID, 'Time\tX1\tX2\tY\n');
fprintf(fileID, '%f\t%f\t%f\t%f\n', data');
fclose(fileID);

disp('Results saved to results_1.txt');
