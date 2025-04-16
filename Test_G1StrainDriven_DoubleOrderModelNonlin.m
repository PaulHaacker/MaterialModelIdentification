% Test Script for G1StrainDriven_DoubleOrderModelNonlin
% Define time vector
clear
close all
T = 100;                % total simulation time [s]
N = 10000;              % number of time steps
time = linspace(0, T, N);
dt = time(2) - time(1);
% Define strain input (e.g., sinusoidal)
strain = sin(2*pi*0.5*time);  % 0.5 Hz, 1% amplitude
% Define model parameters
alpha_1 = .1;
alpha_2 = 0.1;
E_0 = 4;
E_1 = 1;
E_2 = 2;
p_1 = 1;
p_2 = 2;
G = 1;
par = [alpha_1; alpha_2; E_0; E_1; E_2; p_1; p_2; G];
% Initial stress
stress_0 = 0;
% Call the model
stress = G1StrainDriven_DoubleOrderModelNonlin(par, strain, time, stress_0);
% Plotting
figure;
plot(time, strain, '--', 'LineWidth', 1.5); hold on;
plot(time, stress, 'r', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Strain / Stress');
legend('Strain Input', 'Stress Output');
title('Double Order Nonlinear Viscoelastic Response');
grid on;


% First plot data
x = strain(end-1000:end);
y = stress(end-1000:end);

% Compute transformation
a = 6.586; b = 4;
A_mat = [1 a 0 0; 0 0 1 a; -a 1 0 0; 0 0 -a 1];
b_mat = [1; a; -b*a; b];
help2 = A_mat \ b_mat;
A_transf = reshape(help2, 2, 2);
help = [x; y];
transformed_coords = A_transf * help;

% Create a single figure window with two subplots side by side
figure

% Subplot 1: Original strain-stress plot
subplot(1,2,1)
plot(x, y)
xlabel('strain $\varepsilon$', 'Interpreter', 'latex')
ylabel('stress $\sigma$', 'Interpreter', 'latex')
title('Strain-Stress Plot of nonlinear DOM')
set(gca, 'FontSize', 14);

% Subplot 2: Transformed coordinates
subplot(1,2,2)
plot(transformed_coords(1,:), transformed_coords(2,:))
xlabel('Transformed $X$')
ylabel('Transformed $Y$')
title('exaggerated view: linearly transformed data ')
hold on
plot([0 1],[0,a])
help_3 = xlim;
help_3 = help_3(1);
plot([0 help_3], [0 -help_3/a])
legend('$(X~Y)=A(\varepsilon~\sigma)$','eigenvector of $A$ to eigenvalue 1',['eigenvector of $A$ to eigenvalue ',num2str(b)'],'Location','southeast')
set(gca, 'FontSize', 14);
