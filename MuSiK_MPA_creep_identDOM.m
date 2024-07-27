%% Identification of fractional material model from data
clear
close all

load('creep_processed.mat')

% sample the data

number_sample = 1000;
dataStruct = creepPK2_processed;

time = linspace(dataStruct.time(1),dataStruct.time(end),number_sample);
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

alpha_1 = .2;
E_0 = 2500;
E_1 = 500;
p_1 = 50;

alpha_2 = .2;
E_2 = 500;
p_2 = 50;
par_0 =   [ alpha_1;
        alpha_2 ;
        E_0 ;
        E_1 ;
        E_2 ;
        p_1 ;
        p_2 ]; % parameters

[par_lsqnonlin,res] = identify_DoubleOrderModel_creep(...
    time, stress_data,...
    strain_data, par_0);

disp(['Identified parameters: (\\alpha_1,\\alpha_2,E_0,E_1,p_1,E_2,p_2) ='])
disp(num2str((par_lsqnonlin)'))
disp(['residual =', num2str(res)])


% plot results
[t_log,strain_data_log] = samplelog(time, strain_data);
figure
semilogx(time, strain_data,'.-',t_log,strain_data_log,'o',time,G1StressDriven_DoubleOrderModel(par_lsqnonlin,...
    stress_data,time, strain_data(1)))
xlabel('time $t$')
ylabel('strain $\varepsilon(t)$')
title({'Identification of DOM'; ...
        sprintf('Identified parameters: $(\\alpha_1,\\alpha_2,E_0,E_1,p_1,E_2,p_2) = (%s)$', array2strCommas(par_lsqnonlin))})
legend('exp. data','exp. data logarithmically sampled','identified model','Location','southeast')

function str = array2strCommas(array)
    str = sprintf('%.2f, ', array);
    str = str(1:end-2);
end
