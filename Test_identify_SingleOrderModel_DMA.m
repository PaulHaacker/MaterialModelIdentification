%% Identification of fractional material model from data
clear
close all

% generate "data" directly from the model to test the identification - this
% is essentially just a copy of the test script of the model.
alpha = 1/3;
E_0 = 2;
E_1 = 1;
p_1 = 0.5;
par_test =   [alpha;
            E_0;
            E_1;
            p_1]; % parameters in original form
par_norm_test = par2par_norm(par_test);

omega_data = logspace(-5,5,100); % frequency range
ComplexModulus_data = ComplexMod_SingleOrderModel(par_norm_test,omega_data);
storage_data = real(ComplexModulus_data);
loss_data = imag(ComplexModulus_data);

% add noise
noise_lvl = 0.05;
storage_data = addCustomNoise(storage_data,noise_lvl);
loss_data = addCustomNoise(loss_data,noise_lvl);

% parameters to be identified are the normalized parameters which are an
% input to the model "ComplexModulusFcn" and live in R^4, namely
% par_norm ...      (4-by-1)-array of normalized parameters, where 
%                   alpha = par_norm(1) \in (0,1)
%                   b = par_norm(2) = E_1/p_1 > 0
%                   c = par_norm(3) = E_0 + E_1 > 0
%                   d = par_norm(4) = E_0*E_1/p_1> 0

% nonlinear inequality constraint is defined in "nonlincon_identification"

% initial guess
par_0 = ones(4,1);
par_norm0 = par2par_norm(par_0); % normalized parameters

% weight for loss modulus - needed since magnitude of loss data generally
% is much smaller than storage data.
weight_loss = 1;

par_norm_lsqnonlin = ...
    identify_SingleOrderModel_DMA(omega_data, storage_data, loss_data, ...
    par_norm0, weight_loss);

% verify solution:
disp(['True (normalized) parameters:'])
disp(num2str(par_norm_test'))
disp(['Identified (normalized) parameters:'])
disp(num2str(par_norm_lsqnonlin'))

disp(['norm of difference = ',num2str(norm(par_norm_lsqnonlin-par_norm_test))])

% plotting
ComplexModulus_model = ComplexMod_SingleOrderModel(par_norm_lsqnonlin,omega_data);
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


