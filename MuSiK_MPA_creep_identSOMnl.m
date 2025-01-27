%% Identification of fractional material model from data
% clear
close all

load('creep_processed.mat')

% sample the data

dataStruct = creepPK4_processed;

% number_sample = 1000;
number_sample = length(dataStruct.stress);

% end_indx = length(dataStruct.time);
% end_indx = floor(.04*length(dataStruct.time)); % only include ramp
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

% par_0 = [0.06 343 3000 7759 0]; % result of SOM id fixed time step
% par_0 = [0.09336107    0.03936591      33.08783      386.1137   0]; % result of SOM id growing time step
% par_0 = [0.08803826      9.478454      24.87442      147.5583       0];
par_0 = [1 1 1 1 1];
% par_0 = [.3 1000 1000 1000];

tic
[par_lsqnonlin,res] = identify_SingleOrderModelNonlin_creep(...
    time, stress_data,...
    strain_data, par_0);
time_elapsed = toc

disp(['Identified parameters: (alpha E_0 E_1 p_1 G) ='])
disp(num2str((par_lsqnonlin)'))
disp(['residual =', num2str(res)])


%% plot results
[t_log,strain_data_log] = samplelog(time, strain_data);
figure
% semilogx(time, strain_data,'.-',t_log,strain_data_log,'o',time,G1StressDriven_SingleOrderModel(par_norm_lsqnonlin,...
%     stress_data,time, strain_data(1)))
stress_fcn = @(t)interp1(time,stress_data,t,"linear","extrap");
[tModel,strainModel] = G1StressDriven_SingleOrderModelNonlin_growingStepSize(...
    par_lsqnonlin,stress_fcn,[time(1),time(end)],strain_data(1));

semilogx(dataStruct.time, dataStruct.strain*100,'.-',time, strain_data,'.-',...
    t_log,strain_data_log,'ko',...
    tModel,strainModel )
xlabel('time $t$')
ylabel('strain $\varepsilon(t)$')
title({'Identification of SOM nl'; ...
        sprintf('Identified parameters: $(\\alpha,E_0,E_1,p_1,G) = (%s)$', array2strCommas(par_lsqnonlin)); ...
        sprintf('time: %s', num2str(time_elapsed))})
legend('all exp. data','sampled exp. data','exp. data logarithmically sampled','identified model','Location','southeast')

[tModel_initial,strainModel_initial] = G1StressDriven_SingleOrderModelNonlin_growingStepSize(...
    par_0,stress_fcn,[time(1),time(end)],strain_data(1));
figure
loglog(tModel,abs(strainModel-interp1(time,strain_data,tModel,"linear","extrap"))./abs(strainModel))
hold on
loglog(tModel_initial,abs(strainModel_initial-interp1(time,strain_data,tModel,"linear","extrap"))./abs(strainModel_initial))
ylabel('relative error')
xlabel('time')
legend('linear model','nonlinear model')

% figure
% semilogx(time, abs(strain_data-G1StressDriven_SingleOrderModelNonlin_growingStepSize(...
%     par_lsqnonlin,stress_fcn,[time(1),time(end)],strain_data(1))))
% title('error')
% xlabel('time $t$')

function str = array2strCommas(array)
    str = sprintf('%.2f, ', array);
    str = str(1:end-2);
end
