%% Identification of fractional material model from data
clear
close all

%% read creep data
load('creep_processed.mat')

% sample the data

number_sample = 1000;
dataStruct = creepPK2_processed;

time = linspace(dataStruct.time(1),dataStruct.time(end),number_sample);
stress_data = interp1(dataStruct.time,dataStruct.stress,time);
factor_strain = 1; %1.0848; % 1.5;
strain_data = interp1(dataStruct.time,dataStruct.strain,time)*factor_strain;

% verification plot
% Create a new figure
% figure
% tiledlayout('flow')
% nexttile
% semilogx(creepPK2_processed.time, creepPK2_processed.strain, time, strain_data, 'o')
% title('Strain vs Time')
% xlabel('Time')
% ylabel('Strain')
% nexttile
% semilogx(creepPK2_processed.time, creepPK2_processed.stress, time, stress_data, 'o')
% title('Stress vs Time')
% xlabel('Time')
% ylabel('Stress')

% convert to MPa and percent for normalization
stress_data = stress_data/10^6;
strain_data = strain_data*100;

creep.time = time;
creep.stress = stress_data;
creep.strain = strain_data;

% % add noise
% noise_lvl = 0.05;
% storage_data = addCustomNoise(storage_data,noise_lvl);
% loss_data = addCustomNoise(loss_data,noise_lvl);


%% read DMA data

% read data
load('masterSemiSimpleData20deg.mat')
factor = 100; % convert storage and loss in MPa to MPa/%
omega_data = masterData20deg.freq;
storage_data = masterData20deg.storage/factor;
loss_data = masterData20deg.loss/factor;

[omega_data, storage_data, loss_data] = sortArrays(omega_data, storage_data, loss_data);

weight_loss = 10;

DMA.omega = omega_data;
DMA.storage = storage_data;
DMA.loss = loss_data;
DMA.weight_loss = weight_loss;

%% run Simultaneous Identification

% 0.1512739      25.83276      928.5537      4.453071


alpha_1 = .2;
E_0 = 20;
E_1 = 10;
p_1 = 10;

alpha_2 = .05;
E_2 = 100;
p_2 = 50;

par_0 =   [ alpha_1;
        alpha_2 ;
        E_0 ;
        E_1 ;
        E_2 ;
        p_1 ;
        p_2 ]; % parameters

[par_lsqnonlin_frac,res] = identify_DoubleOrderModel_simultaneous(creep, DMA, par_0);

disp(['Identified parameters: (\\alpha_1,\\alpha_2,E_0,E_1,p_1,E_2,p_2) ='])
disp(num2str(par_lsqnonlin_frac'))
disp(['residual =', num2str(res)])


%% plot results
[t_log,strain_data_log] = samplelog(time, strain_data);
figure
semilogx(time, strain_data,'.-',t_log,strain_data_log,'o',time,G1StressDriven_DoubleOrderModel(par_lsqnonlin_frac,...
    stress_data,time, strain_data(1)))
xlabel('time $t$')
ylabel('strain $\varepsilon(t)$')
title({'Identification of DOM'; ...
        sprintf('Identified parameters: $(\\alpha_1,\\alpha_2,E_0,E_1,p_1,E_2,p_2) = (%s)$', array2strCommas(par_lsqnonlin_frac))})
legend('exp. data','exp. data logarithmically sampled','identified model','Location','southeast')
set(gca, 'FontSize', 14)

ComplexModulus_model_frac = ComplexMod_DoubleOrderModel(par_lsqnonlin_frac, omega_data);
storage_model_frac = real(ComplexModulus_model_frac);
loss_model_frac = imag(ComplexModulus_model_frac);

% % fit the whole-order model
% par_norm_lsqnonlin_whole = identify_DoubleOrderModel_DMA(omega_data, storage_data, loss_data, par_norm0, weight_loss);
% par_lsqnonlin_whole = par_norm2par(par_norm_lsqnonlin_whole);
% 
% ComplexModulus_model_whole = ComplexMod_DoubleOrderModel(par_norm_lsqnonlin_whole, omega_data);
% storage_model_whole = real(ComplexModulus_model_whole);
% loss_model_whole = imag(ComplexModulus_model_whole);

% Plot storage modulus
figure;
subplot(2, 1, 1);
hold on;
plot(omega_data, storage_model_frac, '-', 'DisplayName', sprintf('fractional-order Model, Weight Loss $ W_L = %d$, Parameters: $(\\alpha_1,\\alpha_2,E_0,E_1,p_1,E_2,p_2) = (%s)$', weight_loss, array2strCommas(par_lsqnonlin_frac)));
% plot(omega_data, storage_model_whole, '-', 'DisplayName', sprintf('whole-order Model, Weight Loss $ W_L = %d$, Parameters: $(\\alpha,E_0,E_1,p_1) = (%s)$', weight_loss, array2strCommas(par_lsqnonlin_whole)));
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
plot(omega_data, loss_model_frac, '-', 'DisplayName', 'fractional-order Model');
% plot(omega_data, loss_model_whole, '-', 'DisplayName', sprintf('whole-order Model, Weight Loss $ W_L = %d$, Parameters: $(\\alpha,E_0,E_1,p_1) = (%s)$', weight_loss, array2strCommas(par_lsqnonlin_whole)));
plot(omega_data, loss_data, 'o-', 'DisplayName', 'Data');
set(gca,'xscale','log')
set(gca, 'FontSize', 14)
grid on
ylabel('Loss Modulus in MPa')
xlabel('Frequency in Hz')
% legend('show', 'Location','southwest');

% Title for the whole subplot
sgtitle({'Identification of Material Model Parameters using lsqnonlin,';...
    ['i.e. find parameter $p$: $\min_p \sum_{i = data points} ' ...
    '({}^{data}E''_i-{}^{model}E''(p,\omega_i))^2+W_L({}^{data}E''''_i-{}^{model}E''''(p,\omega_i))^2$'];...
    'Data from Mastercurve at $20{}^\circ \mathrm C$ of Temperatures $-10{}^\circ \mathrm C,\dots,+40{}^\circ \mathrm C$'; ...
    ['Model: DOM']});

set(gcf, 'WindowState', 'maximized');

function str = array2strCommas(array)
    str = sprintf('%.2f, ', array);
    str = str(1:end-2);
end

function [sorted_a, sorted_b, sorted_c] = sortArrays(a, b, c)
    combined = [a, b, c];
    sorted_combined = sortrows(combined, 1);

    sorted_a = sorted_combined(:, 1);
    sorted_b = sorted_combined(:, 2);
    sorted_c = sorted_combined(:, 3);
end

%% run DMA id
% % initial guess
% par_0 = ones(4,1);
% par_norm0 = par2par_norm(par_0); % normalized parameters
% 
% % weight for loss modulus - needed since magnitude of loss data generally
% % is much smaller than storage data.
% weight_loss = 1;
% 
% par_norm_lsqnonlin = ...
%     identify_DoubleOrderModel_DMA(omega_data, storage_data, loss_data, ...
%     par_norm0, weight_loss);
% 
% % verify solution:
% disp(['True (normalized) parameters:'])
% disp(num2str(par_norm_test'))
% disp(['Identified (normalized) parameters:'])
% disp(num2str(par_norm_lsqnonlin'))
% 
% disp(['norm of difference = ',num2str(norm(par_norm_lsqnonlin-par_norm_test))])
% 
% % plotting
% ComplexModulus_model = ComplexMod_DoubleOrderModel(par_norm_lsqnonlin,omega_data);
% storage_model = real(ComplexModulus_model);
% loss_model = imag(ComplexModulus_model);
% 
% figure
% tiledlayout('flow')
% nexttile
% plot(omega_data,storage_data,'o-',omega_data,storage_model);
% set(gca,'xscale','log')
% grid on
% legend('data','model')
% ylabel('storage modulus $E"= \Re\{E^\ast\}$')
% xlabel('frequency in Hz')
% nexttile
% plot(omega_data,loss_data,'o-',omega_data,loss_model);
% set(gca,'xscale','log')
% grid on
% legend('data','model')
% ylabel('loss modulus $E""= \Im\{E^\ast\}$')
% xlabel('frequency in Hz')
% 
% 
% function noisy_data = addCustomNoise(data, noise_levels)
%     % Calculate the standard deviation of the data
%     std_dev = std(data);
% 
%     % Calculate the desired noise power
%     noise_power = (noise_levels * std_dev)^2;
% 
%     % Generate random noise with the same length as the data
%     noise = randn(size(data));
% 
%     % Scale the noise to have the desired power
%     scaled_noise = sqrt(noise_power) * noise;
% 
%     % Add the scaled noise to the data
%     noisy_data = data + scaled_noise;
% end


%% run creep id
% % initial guess
% % par_0 = ones(4,1);
% % par_0 = [0.0084621   0.0015443      2.2051   0.0033959]';
% % par_0 = par*100;
% % 
% par_0 = .1*par;
% 
% 
% par_norm0 = par2par_norm(par_0); % normalized parameters
% 
% 
% % par_norm0 = [1,10^3,10^3,10^3]';
% % par_norm0 = [1,1,1,1]';
% 
% [par_norm_lsqnonlin,res] = identify_DoubleOrderModel_creep(time, stress_data,...
%     strain_data, par_norm0);
% 
% % verify solution:
% disp(['True (normalized) parameters:'])
% disp(num2str(par2par_norm(par)'))
% disp(['Identified (normalized) parameters:'])
% disp(num2str(par_norm_lsqnonlin))
% 
% disp(['norm of difference = ',num2str(norm(par_norm_lsqnonlin-par2par_norm(par)))])
% disp(['residual = ',num2str(res)])
% 
% % plot results
% figure
% plot(time,strain_data,'o',time,G1StressDriven_DoubleOrderModel(par_norm_lsqnonlin,...
%     stress_data,time, strain_data(1)))
% xlabel('time $t$')
% ylabel('strain $\varepsilon(t)$')
% title('G1-algorithm applied to $D^\alpha \sigma + b\sigma = cD^\alpha \varepsilon + d\varepsilon$')
% legend('data','identified model')