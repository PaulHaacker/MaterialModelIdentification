%% Script identifying parameters of the 3parameter, fractional, linear material model from measurement data
% data from DMA (dynamic mechanical analysis),
% experiments conducted by Ondrej Farkas and Alexander Lion (Uni-BW Munich)
% material: synthetic resin

clear
close all

% % read data
% load('masterSemiSimpleData20deg.mat')
% omega_data = masterData20deg.freq;
% storage_data = masterData20deg.storage;
% loss_data = masterData20deg.loss;

% read data
load('masterSimpleData20deg.mat')
omega_data = masterSimpleData20deg.freq;
storage_data = masterSimpleData20deg.storage;
loss_data = masterSimpleData20deg.loss;

[omega_data, storage_data, loss_data] = sortArrays(omega_data, storage_data, loss_data);

%% identification and plotting - using different weights

lb = zeros(4,1);
ub = [1 inf inf inf]';

alpha = .2;
E_0 = 2500;
E_1 = 500;
p_1 = 50;
par_0 = [alpha;E_0;E_1;p_1];
par_norm0 = par2par_norm(par_0);

% weight_loss_values = [1, 10, 100]; % Different weight values for loss modulus
weight_loss_values = [1, 5, 10]; % Different weight values for loss modulus

figure;

% Plot storage modulus
subplot(2, 1, 1);
hold on;
for ii = 1:length(weight_loss_values)
    weight_loss = weight_loss_values(ii);

    par_norm_lsqnonlin = identifyMaterialModel(omega_data, storage_data, loss_data, par_norm0, lb, ub, weight_loss, @(par,om)ComplexModulusFcn_3parLinear(par,om), @(p) nonlincon_3parLinear(p));
    par_lsqnonlin = par_norm2par(par_norm_lsqnonlin);

    ComplexModulus_model = ComplexModulusFcn_3parLinear(par_norm_lsqnonlin, omega_data);
    storage_model = real(ComplexModulus_model);

    % Plot storage modulus
    plot(omega_data, storage_model, '-', 'DisplayName', sprintf('Model, Weight Loss $ W_L = %d$, Parameters: $(\\alpha,E_0,E_1,p_1) = (%s)$', weight_loss, array2strCommas(par_lsqnonlin)));
end
plot(omega_data, storage_data, 'o-', 'DisplayName', 'Data');
set(gca,'xscale','log')
set(gca, 'FontSize', 14)
grid on
ylabel('Storage Modulus')
xlabel('Frequency in Hz')
legend('show', 'Location','northwest');

% Plot loss modulus
subplot(2, 1, 2);
hold on;
for ii = 1:length(weight_loss_values)
    weight_loss = weight_loss_values(ii);

    par_norm_lsqnonlin = identifyMaterialModel(omega_data, storage_data, loss_data, par_norm0, lb, ub, weight_loss, @(par,om)ComplexModulusFcn_3parLinear(par,om), @(p) nonlincon_3parLinear(p));
    par_lsqnonlin = par_norm2par(par_norm_lsqnonlin);

    ComplexModulus_model = ComplexModulusFcn_3parLinear(par_norm_lsqnonlin, omega_data);
    loss_model = imag(ComplexModulus_model);

    % Plot loss modulus
    plot(omega_data, loss_model, '-', 'DisplayName', sprintf('Model, Weight Loss $ W_L = %d$, Parameters: $(\\alpha,E_0,E_1,p_1) = (%s)$', weight_loss, array2strCommas(par_lsqnonlin)));
end
plot(omega_data, loss_data, 'o-', 'DisplayName', 'Data');
set(gca,'xscale','log')
set(gca, 'FontSize', 14)
grid on
ylabel('Loss Modulus')
xlabel('Frequency in Hz')
legend('show', 'Location','southwest');

% Title for the whole subplot
sgtitle({'Identification of Material Model Parameters using lsqnonlin,';...
    ['i.e. find parameter $p$: $\min_p \sum_{i = data points} ' ...
    '({}^{data}E''_i-{}^{model}E''(p,\omega_i))^2+W_L({}^{data}E''''_i-{}^{model}E''''(p,\omega_i))^2$'];...
    'Data from Mastercurve at $20{}^\circ \mathrm C$ of Temperatures $-10{}^\circ \mathrm C,\dots,+40{}^\circ \mathrm C$'; ...
    ['Model: $\left(\frac{1}{E_1}{}^\mathrm C D^\alpha_{-\infty}+\frac{1}{p_1}\right)\sigma = ...' ...
    'E_0\left(\left(\frac{1}{E_0}+\frac{1}{E_1}\right){}^\mathrm C D^\alpha_{-\infty}+\frac{1}{p_1}\right)\varepsilon$']});

set(gcf, 'WindowState', 'maximized');

%% identification and grid search over weights
% usage: try out a lot of weights and see whether they converge, i.e.
% compare all identified parameters to the 'last' one, the one using the
% largest weight. This should confirm whether a large weight is good at
% all. Possible Problems: When the weight is chosen astronomically large,
% the identification probably wont give satisfactory results - one expects
% a 'sweetspot'.

lb = zeros(4,1);
ub = [1 inf inf inf]';

alpha = .2;
E_0 = 2500;
E_1 = 500;
p_1 = 50;
par_0 = [alpha;E_0;E_1;p_1];
par_norm0 = par2par_norm(par_0);

% weight_loss_values = 0:.1:50; % Different weight values for loss modulus
weight_loss_values = 1 :20; % Different weight values for loss modulus

par_norm_lsqnonlin = zeros(length(par_norm0),length(weight_loss_values));

for ii = 1:length(weight_loss_values)
    weight_loss = weight_loss_values(ii);

    par_norm_lsqnonlin(:,ii) = identifyMaterialModel(omega_data, storage_data, loss_data, par_norm0, lb, ub, weight_loss, @(par,om)ComplexModulusFcn_3parLinear(par,om), @(p) nonlincon_3parLinear(p));
end

% Compute norm differences
norm_diff_array = vecnorm(par_norm_lsqnonlin-par_norm_lsqnonlin(:,end),2,1);
figure
plot(weight_loss_values,log(norm_diff_array),'o-')

set(gca, 'FontSize', 14)
grid on
ylabel('log of norm of difference of parameters')
xlabel('loss weight $W_L$')

set(gcf, 'WindowState', 'maximized');

%% Comparison with whole-order 3parameter model
% in order to see the benefit of the fractional-order material model, run
% an identification of a whole-order material model by simply constraining
% alpha = 1. As a (dirty) workaround, use a new nonlinear constraint function to
% implement the constraint alpha = 1. A more clean way would be to use the
% linear constraints of lsqnonlin. An even more clean way would be to write
% a new MaterialModelFunction with only 3 parameters and alpha = 1 hard
% coded.

lb = zeros(4,1);
ub = [1 inf inf inf]';

alpha = .2;
E_0 = 2500;
E_1 = 500;
p_1 = 50;
par_0 = [alpha;E_0;E_1;p_1];
par_norm0 = par2par_norm(par_0);

weight_loss = 10; % weight value for loss modulus

% fit the fractional model
par_norm_lsqnonlin_frac = identifyMaterialModel(omega_data, storage_data, loss_data, par_norm0, lb, ub, weight_loss, @(par,om)ComplexModulusFcn_3parLinear(par,om), @(p) nonlincon_3parLinear(p));
par_lsqnonlin_frac = par_norm2par(par_norm_lsqnonlin_frac);

ComplexModulus_model_frac = ComplexModulusFcn_3parLinear(par_norm_lsqnonlin_frac, omega_data);
storage_model_frac = real(ComplexModulus_model_frac);
loss_model_frac = imag(ComplexModulus_model_frac);

% fit the whole-order model
par_norm_lsqnonlin_whole = identifyMaterialModel(omega_data, storage_data, loss_data, par_norm0, lb, ub, weight_loss, @(par,om)ComplexModulusFcn_3parLinear(par,om), @(p) nonlincon_3parLinear_alpha1(p));
par_lsqnonlin_whole = par_norm2par(par_norm_lsqnonlin_whole);

ComplexModulus_model_whole = ComplexModulusFcn_3parLinear(par_norm_lsqnonlin_whole, omega_data);
storage_model_whole = real(ComplexModulus_model_whole);
loss_model_whole = imag(ComplexModulus_model_whole);

% Plot storage modulus
figure;
subplot(2, 1, 1);
hold on;
plot(omega_data, storage_model_frac, '-', 'DisplayName', sprintf('fractional-order Model, Weight Loss $ W_L = %d$, Parameters: $(\\alpha,E_0,E_1,p_1) = (%s)$', weight_loss, array2strCommas(par_lsqnonlin_frac)));
plot(omega_data, storage_model_whole, '-', 'DisplayName', sprintf('whole-order Model, Weight Loss $ W_L = %d$, Parameters: $(\\alpha,E_0,E_1,p_1) = (%s)$', weight_loss, array2strCommas(par_lsqnonlin_whole)));
plot(omega_data, storage_data, 'o-', 'DisplayName', 'Data');
set(gca,'xscale','log')
set(gca, 'FontSize', 14)
grid on
ylabel('Storage Modulus in MPa')
xlabel('Frequency in Hz')
legend('show', 'Location','northwest');

% Plot loss modulus
subplot(2, 1, 2);
hold on;
plot(omega_data, loss_model_frac, '-', 'DisplayName', sprintf('fractional-order Model, Weight Loss $ W_L = %d$, Parameters: $(\\alpha,E_0,E_1,p_1) = (%s)$', weight_loss, array2strCommas(par_lsqnonlin_frac)));
plot(omega_data, loss_model_whole, '-', 'DisplayName', sprintf('whole-order Model, Weight Loss $ W_L = %d$, Parameters: $(\\alpha,E_0,E_1,p_1) = (%s)$', weight_loss, array2strCommas(par_lsqnonlin_whole)));
plot(omega_data, loss_data, 'o-', 'DisplayName', 'Data');
set(gca,'xscale','log')
set(gca, 'FontSize', 14)
grid on
ylabel('Loss Modulus in MPa')
xlabel('Frequency in Hz')
legend('show', 'Location','southwest');

% Title for the whole subplot
sgtitle({'Identification of Material Model Parameters using lsqnonlin,';...
    ['i.e. find parameter $p$: $\min_p \sum_{i = data points} ' ...
    '({}^{data}E''_i-{}^{model}E''(p,\omega_i))^2+W_L({}^{data}E''''_i-{}^{model}E''''(p,\omega_i))^2$'];...
    'Data from Mastercurve at $20{}^\circ \mathrm C$ of Temperatures $-10{}^\circ \mathrm C,\dots,+40{}^\circ \mathrm C$'; ...
    ['Model: $\left(\frac{1}{E_1}{}^\mathrm C D^\alpha_{-\infty}+\frac{1}{p_1}\right)\sigma = ' ...
    'E_0\left(\left(\frac{1}{E_0}+\frac{1}{E_1}\right){}^\mathrm C D^\alpha_{-\infty}+\frac{1}{p_1}\right)\varepsilon$']});

set(gcf, 'WindowState', 'maximized');


%% functions
function [sorted_a, sorted_b, sorted_c] = sortArrays(a, b, c)
    combined = [a, b, c];
    sorted_combined = sortrows(combined, 1);

    sorted_a = sorted_combined(:, 1);
    sorted_b = sorted_combined(:, 2);
    sorted_c = sorted_combined(:, 3);
end

function str = array2strCommas(array)
    str = sprintf('%.2f, ', array);
    str = str(1:end-2);
end
