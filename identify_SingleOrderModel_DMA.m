function par_norm_lsqnonlin = identify_SingleOrderModel_DMA(omega_data, storage_data, loss_data, initial_guess,weight_loss)
% identifyMaterialModel Identifies parameters of a material model using experimental data.
%
%   par_norm_lsqnonlin = identifyParameters(omega_data, storage_data, loss_data, initial_guess, weight_loss)
%
% Inputs:
%   omega_data           - Frequencies at which the material properties are evaluated.
%   storage_data         - Experimental storage modulus data in the frequency domain.
%   loss_data            - Experimental loss modulus data in the frequency domain.
%   initial_guess        - Initial guess for the (normalized) parameters to be identified.
%   weight_loss          - Weighting factor for the loss modulus data.
%
% Output:
%   par_norm_lsqnonlin   - Identified normalized parameters obtained through optimization.

% The parameters are subject to constraints, namely they are partially
% bounded:
lb = zeros(4,1);
ub = [1 inf inf inf]';


% Define options for lsqnonlin
options = optimoptions("lsqnonlin","Algorithm","interior-point",'Display','off');

% Call lsqnonlin to optimize parameters
par_norm_lsqnonlin = lsqnonlin(@(p) diffFcn(storage_data, loss_data, weight_loss, ComplexMod_SingleOrderModel(p,omega_data)), ...
    initial_guess, lb, ub, [], [], [], [], @(p) nonlincon_SingleOrderModel(p), options);

end

function diff = diffFcn(storage_data,loss_data,weight_loss,ComplexMod)
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
