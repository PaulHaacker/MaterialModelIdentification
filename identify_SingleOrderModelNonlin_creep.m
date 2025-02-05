function [par_norm_lsqnonlin,res] = identify_SingleOrderModelNonlin_creep(time, stress_data, strain_data, initial_guess)
% Identifies parameters of a material model using experimental data.
%
%   par_norm_lsqnonlin = identify_SingleOrderModel_creep(time, stress_data, strain_data, initial_guess)
%
% Inputs:

% Output:
%   par_norm_lsqnonlin   - Identified normalized parameters obtained through optimization.

warning off
% warning('off', 'MATLAB:interp1:NaNinY');

% The parameters are subject to constraints, namely they are partially
% bounded:
lb = zeros(5,1);
lb(5) = -10^4*0;
% lb(1) =0.03;
% lb(1) =0.15;
% ub = [1 inf inf inf]';
ub = [1 ones(1,4)*10^4]';
% ub(1) =  .1;
% ub(1) =  .18;
% ub(2) = 100;

% Define options for lsqnonlin
% options = optimoptions("lsqnonlin","Algorithm","interior-point",'Display','off');
options = optimoptions("lsqnonlin","Algorithm","interior-point",'Display','off','OptimalityTolerance',1e-9, ...
    'StepTolerance',0,'ConstraintTolerance',0);

% Call lsqnonlin to optimize parameters
% [par_norm_lsqnonlin,res] = lsqnonlin(@(p) GrowingTimeStepDiff(p,strain_data,stress_data,time), ...
%     initial_guess, lb, ub, [], [], [], [], @(p) nonlincon_SingleOrderModel(p), options);
[par_norm_lsqnonlin,res] = lsqnonlin(@(p) GrowingTimeStepDiff(p,strain_data,stress_data,time), ...
    initial_guess, lb, ub,  options);
% [par_norm_lsqnonlin,res] = lsqnonlin(@(p) LogTimeStepDiff(p,strain_data,stress_data,time), ...
%     initial_guess, lb, ub, [], [], [], [], @(p) nonlincon_SingleOrderModel(p), options);
% [par_norm_lsqnonlin,res] = lsqnonlin(@(p) LogTimeStepDiff(p,strain_data,stress_data,time), ...
%     initial_guess, lb, ub, [], [], [], [], @(p) nonlincon_SingleOrderModel_alpha1(p), options);
% [par_norm_lsqnonlin,res] = lsqnonlin(@(p) EquidistTimeStepDiff(p,strain_data,stress_data,time), ...
%     initial_guess, lb, ub, [], [], [], [], @(p) nonlincon_SingleOrderModel(p), options);
% [par_norm_lsqnonlin,res] = lsqnonlin(@(p) EquidistTimeStepDiff(p,strain_data,stress_data,time), ...
%     initial_guess, lb, ub);
initial_res = sum(GrowingTimeStepDiff(initial_guess,strain_data,stress_data,time).^2);
if initial_res < res
    warning('residual has increased - local minimum found.')
end
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
    % Debugging: Check if diff_logdist contains NaN or Inf
    if any(isnan(diff_logdist)) || any(isinf(diff_logdist))
        error('diff_logdist contains NaN or Inf values');
    end
end

function diff_GrowDist = GrowingTimeStepDiff(p,strain_data,stress_data,time)
    stress_fcn = @(t) interp1(time,stress_data, t,'linear','extrap');
    % compute gruenwald time stepping
    [t_vec_sim, strain_vec_sim] = G1StressDriven_SingleOrderModelNonlin_growingStepSize(p,stress_fcn,[time(1),time(end)],strain_data(1));
    % diff_GrowDist = sqrt(abs(strain_vec_sim - interp1(time, strain_data,t_vec_sim,'linear','extrap'))./abs(strain_vec_sim)); % relative error
    diff_GrowDist = (strain_vec_sim - interp1(time, strain_data,t_vec_sim,'linear','extrap')); % absolute error
end