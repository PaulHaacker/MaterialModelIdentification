%% Identification of fractional material model from data
clear
close all

load('creep_processed.mat')

% sample the data

dataStruct = creepPK4_processed;

% number_sample = 1000;
number_sample = length(dataStruct.stress);

% end_indx = length(dataStruct.time);
end_indx = floor(.8*length(dataStruct.time)); % cut off last bit of data which is relaxation

% time = linspace(dataStruct.time(1),dataStruct.time(end),number_sample);
time = linspace(dataStruct.time(1),dataStruct.time(end_indx),number_sample);
time = time; % modify for logarithmic time scale
stress_data = interp1(dataStruct.time,dataStruct.stress,time);
strain_data = interp1(dataStruct.time,dataStruct.strain,time);

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
% par_0 = [.3 1000 1000 1000];
% par_0 = [1 1 1 1];
par_0 = [0.06 343 3000 7759];
% par_0 = [.3 1000 1000 1000];
% par_0 = par*1.1;
par_norm0 = par2par_norm(par_0); % normalized parameters

[par_norm_lsqnonlin,res] = identify_SingleOrderModel_creep(...
    time, stress_data,...
    strain_data, par_norm0);

disp(['Identified parameters: (alpha E_0 E_1 p_1) ='])
disp(num2str(par_norm2par(par_norm_lsqnonlin)'))
disp(['residual =', num2str(res)])


% plot results
[t_log,strain_data_log] = samplelog(time, strain_data);
figure
% semilogx(time, strain_data,'.-',t_log,strain_data_log,'o',time,G1StressDriven_SingleOrderModel(par_norm_lsqnonlin,...
%     stress_data,time, strain_data(1)))
semilogx(dataStruct.time, dataStruct.strain*100,'.-',time, strain_data,'.-',t_log,strain_data_log,'ko',time,G1StressDriven_SingleOrderModel(par_norm_lsqnonlin,...
    stress_data,time, strain_data(1)))
xlabel('time $t$')
ylabel('strain $\varepsilon(t)$')
title({'Identification of $D^\alpha \sigma + b\sigma = cD^\alpha \varepsilon + d\varepsilon$'; ...
        sprintf('Identified parameters: $(\\alpha,b,c,d) = (%s)$', array2strCommas(par_norm_lsqnonlin))})
legend('all exp. data','sampled exp. data','exp. data logarithmically sampled','identified model','Location','southeast')

figure
semilogx(time, abs(strain_data-G1StressDriven_SingleOrderModel(par_norm_lsqnonlin,...
    stress_data,time, strain_data(1))))
title('error')
xlabel('time $t$')
function str = array2strCommas(array)
    str = sprintf('%.2f, ', array);
    str = str(1:end-2);
end
