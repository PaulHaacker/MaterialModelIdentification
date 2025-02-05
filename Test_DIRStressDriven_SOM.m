% Test script for DIRStressDriven_SOM

% Define parameters
par = [1/3 2 1 1/2];  % [alpha, E0, E1, p1]

% Define stress as a constant function of time
stress = @(t) 1;  % Constant stress signal

% Define time span and initial strain
tspan = [0, 1];  % From t = 0 to t = 5
strain_0 = 0;    % Initial strain condition

% Call DIRStressDriven_SOM
[t_vec, strain_vec] = DIRStressDriven_SOM(par, stress, tspan, strain_0);

% Plot the results
figure;
plot(t_vec, strain_vec, 'b-', 'LineWidth', 1.5);
xlabel('Time (t)', 'FontWeight', 'bold');
ylabel('Strain (\epsilon(t))', 'FontWeight', 'bold');
title('Strain Response in Single-Order Viscoelastic Model', 'FontWeight', 'bold');
grid on;
box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
