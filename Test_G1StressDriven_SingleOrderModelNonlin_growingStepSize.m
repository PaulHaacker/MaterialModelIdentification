clear
close all

alpha = 1/3;
E_0 = 2;
E_1 = 1;
p_1 = 0.5;

% Define G values to iterate over
G_values = [-1, 0, 1];

tspan = [0,1]; % time span
stress = @(t) t.^0;
strain_0 = 0;

% Prepare figure for plotting
figure
hold on % allows multiple plots on the same figure

% Loop over different G values
for G = G_values
    % Update parameters with current G
    par =   [alpha;
            E_0;
            E_1;
            p_1;
            G]; % parameters
    
    % Call the function with updated parameters
    [t_vec, strain_vec] = G1StressDriven_SingleOrderModelNonlin_growingStepSize(par, stress, tspan, strain_0);
    
    % Plot the result for this G value
    plot(t_vec, strain_vec, 'o-', 'DisplayName', ['G = ', num2str(G)])
end

% Label the plot
xlabel('time $t$', 'Interpreter', 'latex')
ylabel('strain $\varepsilon(t)$', 'Interpreter', 'latex')
title('G1-algorithm with varying $G$', 'Interpreter', 'latex')

% Add a legend to distinguish the different G values
legend('Location','southeast')

% Finalize plot settings
hold off
set(gca, 'FontSize', 16);
