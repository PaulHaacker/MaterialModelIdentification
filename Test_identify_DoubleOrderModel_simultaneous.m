%% Identification of fractional material model from data
clear
close all

%% generate creep data
% generate "data" directly from the model to test the identification - this
% is essentially just a copy of the test script of the model.
alpha_1 = .2;
E_0 = 2500;
E_1 = 500;
p_1 = 50;

alpha_2 = .2;
E_2 = 500;
p_2 = 50;
par_test =   [ alpha_1;
        alpha_2 ;
        E_0 ;
        E_1 ;
        E_2 ;
        p_1 ;
        p_2 ]; % parameters

time = 0:.01:10; % time array
% stress_data = sin(time.^2); % 
% stress_data = min(time / 0.5, 1); % ramp up to 1 at 0.5 seconds and stay constant
stress_data = ones(size(time)); % ramp up to 1 at 0.5 seconds and stay constant

strain_0 = 0;

strain_data = G1StressDriven_DoubleOrderModel((par_test),stress_data,time,strain_0);

creep.time = time;
creep.stress = stress_data;
creep.strain = strain_data;

% % add noise
% noise_lvl = 0.05;
% storage_data = addCustomNoise(storage_data,noise_lvl);
% loss_data = addCustomNoise(loss_data,noise_lvl);


%% generate DMA data
% generate "data" directly from the model to test the identification - this
% is essentially just a copy of the test script of the model.

omega_data = logspace(-5,5,100); % frequency range
ComplexModulus_data = ComplexMod_DoubleOrderModel(par_test,omega_data);
storage_data = real(ComplexModulus_data);
loss_data = imag(ComplexModulus_data);

weight_loss = 10;

% % add noise
% noise_lvl = 0.05;
% storage_data = addCustomNoise(storage_data,noise_lvl);
% loss_data = addCustomNoise(loss_data,noise_lvl);

DMA.omega = omega_data;
DMA.storage = storage_data;
DMA.loss = loss_data;
DMA.weight_loss = weight_loss;

%% run Simultaneous Identification
par_0 = ones(7,1);
[par_lsqnonlin,res] = identify_DoubleOrderModel_simultaneous(creep, DMA, par_0);

% verify solution:
disp(['True parameters:'])
disp(num2str((par_test)'))
disp(['Identified parameters:'])
disp(num2str(par_lsqnonlin'))

disp(['norm of difference = ',num2str(norm(par_lsqnonlin-(par_test)))])
disp(['residual = ',num2str(res)])

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