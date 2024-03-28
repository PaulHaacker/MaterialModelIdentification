function par_norm_lsqnonlin = identifyMaterialModel(omega_data, storage_data, loss_data, initial_guess, lb, ub, weight_loss, ComplexModulusFcn, nonlincon)
% identifyMaterialModel Identifies parameters of a material model using experimental data.
%
%   par_norm_lsqnonlin = identifyParameters(omega_data, storage_data, loss_data, initial_guess, lb, ub, weight_loss, ComplexModulusFcn)
%
% Inputs:
%   omega_data           - Frequencies at which the material properties are evaluated.
%   storage_data         - Experimental storage modulus data in the frequency domain.
%   loss_data            - Experimental loss modulus data in the frequency domain.
%   initial_guess        - Initial guess for the (normalized) parameters to be identified.
%   lb                   - Lower bounds for the parameters.
%   ub                   - Upper bounds for the parameters.
%   weight_loss          - Weighting factor for the loss modulus data.
%   ComplexModulusFcn    - Function handle representing a model of the complex modulus. 
%                          This function should take normalized parameters and frequencies
%                          as inputs and return the complex modulus.
%
% Output:
%   par_norm_lsqnonlin   - Identified normalized parameters obtained through optimization.

    % Define options for lsqnonlin
    options = optimoptions("lsqnonlin","Algorithm","interior-point",'Display','off');

    % Call lsqnonlin to optimize parameters
    par_norm_lsqnonlin = lsqnonlin(@(p) diffFcn(p, omega_data, storage_data, loss_data, weight_loss, ComplexModulusFcn), ...
        initial_guess, lb, ub, [], [], [], [], @(p) nonlincon(p), options);
end

function diff = diffFcn(par_norm,omega,storage_data,loss_data,weight_loss,ComplexModulusFcn)
% difference of material model in frequency domain and data. for
% identification pruposes, namely to be handed to lsqnonlin

% Inputs:
%   par_norm            - Normalized parameters of the material model
%   omega               - frequencies at which the material properties are evaluated
%   storage_data        - Experimental storage modulus data in frequency domain
%   loss_data           - Experimental loss modulus data in frequency domain
%   weight_loss         - Weighting factor for the loss modulus data
%   ComplexModulusFcn   - function of par_norm and omega, representing a
%                         model of the complex modulus
%
% Outputs:
%   diff          - Difference between the experimental data and the material model
%                  It is a column vector with two components:
%                  1. Difference between experimental and model storage modulus data
%                  2. Weighted difference between experimental and model loss modulus data

ComplModulusModel = ComplexModulusFcn(par_norm,omega);
storage_model = real(ComplModulusModel);
loss_model = imag(ComplModulusModel);

diff = [storage_data - storage_model;
        weight_loss*(loss_data - loss_model)];
end
