clear
close all
% Define parameters
par = [0.3, 1000, 500, 10, 200]; % Example parameters

% Define time span
tspan = [0, 10]; % 10-second simulation

% Define initial strain
strain_0 = 0;

% Define time vector
time_steps = linspace(tspan(1), tspan(2), 3000); % Equidistant time points

% Define sinusoidal stress input
stress_amplitude = 100;
stress_frequency = 1; % 1 Hz
stress = @(t) stress_amplitude * sin(2 * pi * stress_frequency * t);

% Call the function
[t_vec, strain_vec] = G1StressDriven_SingleOrderModelNonlin(par, stress, time_steps, strain_0);

% Plot results
figure;
subplot(2,1,1);
plot(t_vec, stress(t_vec), 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Stress (Pa)');
title('Sinusoidal Stress Input');
grid on;

subplot(2,1,2);
plot(t_vec, strain_vec, 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Strain');
title('Resulting Strain Response');
grid on;
