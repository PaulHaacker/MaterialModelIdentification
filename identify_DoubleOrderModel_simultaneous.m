function [par_lsqnonlin,res] = identify_DoubleOrderModel_simultaneous(creep, DMA, initial_guess)
% Identifies parameters of a material model using experimental data.
%
%   par_norm_lsqnonlin = identify_DoubleOrderModel_creep(time, stress_data, strain_data, initial_guess)
%
% Inputs: creep and DMA data
time = creep.time;
stress_data = creep.stress;
strain_data = creep.strain;

omega_data = DMA.omega;
storage_data = DMA.storage;
loss_data = DMA.loss;
weight_loss = DMA.weight_loss;

% Output:
%   par_norm_lsqnonlin   - Identified normalized parameters obtained through optimization.

warning off
% warning('off', 'MATLAB:interp1:NaNinY');

% The parameters are subject to constraints
lb = zeros(7,1);
% lb(1:2) = .02;
% lb(1) =0.15;
% ub = [1 inf inf inf]';
ub = ones(7,1);
% ub(1:2) = .2;
ub(3:end) = ub(3:end)*10^4;
% ub = [1 ones(1,3)*10^4]';
% ub(1) =  .18;
% ub(2) = 100;

% Define options for lsqnonlin
options = optimoptions("lsqnonlin","Algorithm","interior-point",'Display','off');
% options = optimoptions("lsqnonlin","Algorithm","interior-point",'Display','off','OptimalityTolerance',1e-9, ...
%     'StepTolerance',0,'ConstraintTolerance',0);

% Optimize creep alone
CreepDiff = @(p) LogTimeStepDiff(p,strain_data,stress_data,time);
[par_creep,res_creep] = lsqnonlin(CreepDiff, ...
    initial_guess, lb, ub, options);
% [par_norm_lsqnonlin,res] = lsqnonlin(@(p) LogTimeStepDiff(p,strain_data,stress_data,time), ...
%     initial_guess, lb, ub, [], [], [], [], @(p) nonlincon_DoubleOrderModel_alpha1(p), options);
% [par_norm_lsqnonlin,res] = lsqnonlin(@(p) EquidistTimeStepDiff(p,strain_data,stress_data,time), ...
%     initial_guess, lb, ub, [], [], [], [], @(p) nonlincon_DoubleOrderModel(p), options);
% [par_norm_lsqnonlin,res] = lsqnonlin(@(p) EquidistTimeStepDiff(p,strain_data,stress_data,time), ...
%     initial_guess, lb, ub);

% Optimize DMA alone
DMADiff = @(p) reshape(...
    diffFcnDMA(storage_data, loss_data, weight_loss, ComplexMod_DoubleOrderModel(p,omega_data)),...
    [2*length(omega_data),1]);
[par_DMA,res_DMA] = lsqnonlin(DMADiff, ...
    initial_guess, lb, ub, options);

% weight_creep = res_DMA/res_creep;
% weight_creep = 1.97*10^9;
% weight_creep = 1*10^9;
% weight_creep = 1*10^3;
weight_creep = 1*10^3;
% weight_DMA = 1;
% weight_DMA = res_creep/res_DMA;

[par_lsqnonlin,res] = lsqnonlin(@(p) [DMADiff(p);sqrt(weight_creep)*CreepDiff(p)'], ...
    par_creep, lb, ub, options);
% [par_lsqnonlin,res] = lsqnonlin(@(p) [DMADiff(p);sqrt(weight_creep)*CreepDiff(p)'], ...
%     initial_guess, lb, ub, options);
% [par_lsqnonlin,res] = lsqnonlin(@(p) [sqrt(weight_DMA)*DMADiff(p);CreepDiff(p)'], ...
%     par_DMA, lb, ub, options);
end

function diff = diffFcnDMA(storage_data,loss_data,weight_loss,ComplexMod)
% difference of material model in frequency domain and data. for
% identification pruposes, namely to be handed to lsqnonlin

% Inputs:
%   storage_data        - Experimental storage modulus data in frequency domain
%   loss_data           - Experimental loss modulus data in frequency domain
%   weight_loss         - Weighting factor for the loss modulus data
%   ComplexMod          - function of par_norm and omega, representing a
%                         model of the complex modulus
%
% Outputs:
%   diff          - Difference between the experimental data and the material model
%                  It is a column vector with two components:
%                  1. Difference between experimental and model storage modulus data
%                  2. Weighted difference between experimental and model loss modulus data

storage_model = real(ComplexMod);
loss_model = imag(ComplexMod);

diff = [storage_data - storage_model;
        weight_loss*(loss_data - loss_model)];
end

function diff_equidist = EquidistTimeStepDiff(p,strain_data,stress_data,time)
    % compute gruenwald time stepping
    diff_equidist = strain_data - G1StressDriven_DoubleOrderModel(p,stress_data,time,strain_data(1));
end

function diff_logdist = LogTimeStepDiff(p,strain_data,stress_data,time)
    % compute gruenwald time stepping
    diff_equidist = strain_data - G1StressDriven_DoubleOrderModel(p,stress_data,time,strain_data(1));
    % sample logarithmically
    [~,diff_logdist] = samplelog(time, diff_equidist);
    % Debugging: Check if diff_logdist contains NaN or Inf
    if any(isnan(diff_logdist)) || any(isinf(diff_logdist))
        error('diff_logdist contains NaN or Inf values');
    end
end