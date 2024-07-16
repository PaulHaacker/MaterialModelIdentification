function par_norm_lsqnonlin = identify_SingleOrderModel_creep(time, stress_data, strain_data, initial_guess)
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
ub = [1 inf inf inf]';

% Define options for lsqnonlin
options = optimoptions("lsqnonlin","Algorithm","interior-point",'Display','off');

% Call lsqnonlin to optimize parameters
par_norm_lsqnonlin = lsqnonlin(@(p) strain_data - G1StressDriven_SingleOrderModel(p,stress_data,time,strain_data(1)), ...
    initial_guess, lb, ub, [], [], [], [], @(p) nonlincon_SingleOrderModel(p), options);
end