%% Identification of fractional material model from data
clear
close all

% generate "data" directly from the model to test the identification - this
% is essentially just a copy of the test script of the model.
alpha = 1/3;
E_0 = 2;
E_1 = 1;
p_1 = 0.5;
par =   [alpha;
        E_0;
        E_1;
        p_1]; % parameters

time = 0:.01:10; % time array
% stress_data = sin(time.^2); % 
% stress_data = min(time / 0.5, 1); % ramp up to 1 at 0.5 seconds and stay constant
stress_fcn = @(t) ones(size(t));
stress_data = stress_fcn(time); % step input

strain_0 = 0;

% stress_fcn = @(t) interp1(time,stress_data, t);
% [t_vec_data, strain_data] = G1StressDriven_SingleOrderModel_growingStepSize(par2par_norm(par),stress_fcn,[time(1),time(end)],strain_0);

strain_data = G1StressDriven_SingleOrderModel(par2par_norm(par),stress_data,time,strain_0);

% % add noise
% noise_lvl = 0.05;
% storage_data = addCustomNoise(storage_data,noise_lvl);
% loss_data = addCustomNoise(loss_data,noise_lvl);

% parameters to be identified are the normalized parameters which are an
% input to the model "ComplexModulusFcn" and live in R^4, namely
% par_norm ...      (4-by-1)-array of normalized parameters, where 
%                   alpha = par_norm(1) \in (0,1)
%                   b = par_norm(2) = E_1/p_1 > 0
%                   c = par_norm(3) = E_0 + E_1 > 0
%                   d = par_norm(4) = E_0*E_1/p_1> 0

% initial guess
% par_0 = ones(4,1);
% par_0 = [0.0084621   0.0015443      2.2051   0.0033959]';
% par_0 = par*100;
% 
par_0 = 0.9*par;


par_norm0 = par2par_norm(par_0); % normalized parameters


% par_norm0 = [1,10^3,10^3,10^3]';
% par_norm0 = [1,1,1,1]';

tic
[par_norm_lsqnonlin,res] = identify_SingleOrderModel_creep(time, stress_data,...
    strain_data, par_norm0);
time_elapsed = toc;

% verify solution:
disp(['True (normalized) parameters:'])
disp(num2str(par2par_norm(par)'))
disp(['Identified (normalized) parameters:'])
disp(num2str(par_norm_lsqnonlin))

disp(['norm of difference = ',num2str(norm(par_norm_lsqnonlin-par2par_norm(par)))])
disp(['residual = ',num2str(res)])

% plot results
figure
plot(time,strain_data,'o',time,G1StressDriven_SingleOrderModel(par_norm_lsqnonlin,...
    stress_data,time, strain_data(1)))
xlabel('time $t$')
ylabel('strain $\varepsilon(t)$')
title(['dry run identification results of SOM, time elapsed: ', num2str(time_elapsed)])
legend('data','identified model')


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


