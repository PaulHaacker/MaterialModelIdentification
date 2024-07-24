%% Identification of fractional material model from data
clear
close all

load('creep_processed.mat')

% sample the data

number_sample = 1000;

time = linspace(creepPK2_processed.time(1),creepPK2_processed.time(end),number_sample);
time = time +1; % modify for logarithmic time scale
stress_data = interp1(creepPK2_processed.time,creepPK2_processed.stress,time);
strain_data = interp1(creepPK2_processed.time,creepPK2_processed.strain,time);

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

% parameters to be identified are the normalized parameters which are an
% input to the model "ComplexModulusFcn" and live in R^4, namely
% par_norm ...      (4-by-1)-array of normalized parameters, where 
%                   alpha = par_norm(1) \in (0,1)
%                   b = par_norm(2) = E_1/p_1 > 0
%                   c = par_norm(3) = E_0 + E_1 > 0
%                   d = par_norm(4) = E_0*E_1/p_1> 0

% initial guess
% par_0 = [.16 2222 1000 1000];
par_0 = [.5 1000 1000 1000];
% par_0 = [.3 1000 1000 1000];
% par_0 = par*1.1;
par_norm0 = par2par_norm(par_0); % normalized parameters

par_norm_lsqnonlin = identify_SingleOrderModel_creep(...
    time, stress_data,...
    strain_data, par_norm0);

disp(['Identified parameters: (alpha E_0 E_1 p_1) ='])
disp(num2str(par_norm2par(par_norm_lsqnonlin)'))

% plot results
figure
semilogx(time,strain_data,'o',time,G1StressDriven_SingleOrderModel(par_norm_lsqnonlin,...
    stress_data,time, strain_data(1)))
xlabel('time $t$')
ylabel('strain $\varepsilon(t)$')
title('G1-algorithm applied to $D^\alpha \sigma + b\sigma = cD^\alpha \varepsilon + d\varepsilon$')
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


