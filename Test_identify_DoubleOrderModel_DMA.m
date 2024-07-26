%% Identification of fractional material model from data
clear
close all


%% generate data
alpha_1 = .5;
alpha_2 = .2;
E_0 = 100;
E_1 = 20;
E_2 = 20;
p_1 = 40;
p_2 = 40;
par_test =   [ alpha_1;
        alpha_2 ;
        E_0 ;
        E_1 ;
        E_2 ;
        p_1 ;
        p_2 ]; % parameters

omega_data = logspace(-5,5,100); % frequency range
ComplexModulus_data = ComplexMod_DoubleOrderModel(par_test,omega_data);
storage_data = real(ComplexModulus_data);
loss_data = imag(ComplexModulus_data);

% % add noise
% noise_lvl = 0.05;
% storage_data = addCustomNoise(storage_data,noise_lvl);
% loss_data = addCustomNoise(loss_data,noise_lvl);

%% identification

% initial guess
par_0 = ones(7,1);
% par_0 = par_test;

% weight for loss modulus - needed since magnitude of loss data generally
% is much smaller than storage data.
weight_loss = 1;

par_norm_lsqnonlin = ...
    identify_DoubleOrderModel_DMA(omega_data, storage_data, loss_data, ...
    par_0, weight_loss);

% verify solution:
disp(['True (normalized) parameters:'])
disp(num2str(par_test'))
disp(['Identified (normalized) parameters:'])
disp(num2str(par_norm_lsqnonlin'))

disp(['norm of difference = ',num2str(norm(par_norm_lsqnonlin-par_test))])

% plotting
ComplexModulus_model = ComplexMod_DoubleOrderModel(par_norm_lsqnonlin,omega_data);
storage_model = real(ComplexModulus_model);
loss_model = imag(ComplexModulus_model);

figure
tiledlayout('flow')
nexttile
plot(omega_data,storage_data,'o-',omega_data,storage_model);
set(gca,'xscale','log')
grid on
legend('data','model')
ylabel('storage modulus $E"= \Re\{E^\ast\}$')
xlabel('frequency in Hz')
nexttile
plot(omega_data,loss_data,'o-',omega_data,loss_model);
set(gca,'xscale','log')
grid on
legend('data','model')
ylabel('loss modulus $E""= \Im\{E^\ast\}$')
xlabel('frequency in Hz')


function noisy_data = addCustomNoise(data, noise_levels)
    % Calculate the standard deviation of the data
    std_dev = std(data);

    % Calculate the desired noise power
    noise_power = (noise_levels * std_dev)^2;

    % Generate random noise with the same length as the data
    noise = randn(size(data));

    % Scale the noise to have the desired power
    scaled_noise = sqrt(noise_power) * noise;

    % Add the scaled noise to the data
    noisy_data = data + scaled_noise;
end


