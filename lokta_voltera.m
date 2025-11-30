function yprime = lv(~, y)
    % Lotka-Volterra eqns
    a = 0.5471;
    b = 0.0281;
    r = 0.8439;
    c = 0.0266;
    
    yprime = [a * y(1) - b * y(1) * y(2);
              -r * y(2) + c * y(1) * y(2)];
end

% Initial conditions
y0 = [30; 4];

% Time span (0 to 20 years)
tspan = [0 20];

% Solve 
[t, y] = ode45(@lv, tspan, y0);

% Plot 
figure;
plot(t, y(:,1), 'o');
title('Prey Population Over Time');
xlabel('Time (years)');
ylabel('Prey Population');
grid on;

 