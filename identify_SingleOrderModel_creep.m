function [par_norm_lsqnonlin,res] = identify_SingleOrderModel_creep(time, stress_data, strain_data, initial_guess)
% Identifies parameters of a material model using experimental data.
%
%   par_norm_lsqnonlin = identify_SingleOrderModel_creep(time, stress_data, strain_data, initial_guess)
%
% Inputs:

% Output:
%   par_norm_lsqnonlin   - Identified normalized parameters obtained through optimization.

% The parameters are subject to constraints, namely they are partially
% bounded:
lb = zeros(4,1);
% ub = [1 inf inf inf]';
ub = [1 ones(1,3)*10^4]';

% Define options for lsqnonlin
% options = optimoptions("lsqnonlin","Algorithm","interior-point",'Display','off');
options = optimoptions("lsqnonlin","Algorithm","interior-point",'Display','off','OptimalityTolerance',1e-9, ...
    'StepTolerance',0,'ConstraintTolerance',0);

% Call lsqnonlin to optimize parameters
% [par_norm_lsqnonlin,res] = lsqnonlin(@(p) LogTimeStepDiff(p,strain_data,stress_data,time), ...
%     initial_guess, lb, ub, [], [], [], [], @(p) nonlincon_SingleOrderModel(p), options);
[par_norm_lsqnonlin,res] = lsqnonlin(@(p) EquidistTimeStepDiff(p,strain_data,stress_data,time), ...
    initial_guess, lb, ub, [], [], [], [], @(p) nonlincon_SingleOrderModel(p), options);
% [par_norm_lsqnonlin,res] = lsqnonlin(@(p) EquidistTimeStepDiff(p,strain_data,stress_data,time), ...
%     initial_guess, lb, ub);
end


function diff_equidist = EquidistTimeStepDiff(p,strain_data,stress_data,time)
    % compute gruenwald time stepping
    diff_equidist = strain_data - G1StressDriven_SingleOrderModel(p,stress_data,time,strain_data(1));
end

function diff_logdist = LogTimeStepDiff(p,strain_data,stress_data,time)
    % compute gruenwald time stepping
    diff_equidist = strain_data - G1StressDriven_SingleOrderModel(p,stress_data,time,strain_data(1));
    % sample logarithmically
    [~,diff_logdist] = samplelog(time, diff_equidist);
end