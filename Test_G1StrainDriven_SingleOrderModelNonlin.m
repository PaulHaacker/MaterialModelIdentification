% Test script for G1StrainDriven_SingleOrderModelNonlin
% Define parameters
par = [0.3, 1, 2, 3, 1];  % Example parameters: [alpha, E0, E1, p1, G]
time = linspace(0, 100, 1000);  % Time vector from 0 to 10s with 100 points
strain = @(t)  ones(size(t));  % Example sinusoidal strain function
stress_0 = 0;  % Initial stress condition
% Call the function
[time_out, stress_vec] = G1StrainDriven_SingleOrderModelNonlin(par, strain, time, stress_0);
% Plot results
figure;
plot(time_out, stress_vec, 'b-', 'LineWidth', 2);
hold on;
plot(time, strain(time), 'r--', 'LineWidth', 2);
hold off;
grid on;
legend('Stress', 'Strain');
xlabel('Time (s)');
ylabel('Stress / Strain');
title('Viscoelastic Material Response');