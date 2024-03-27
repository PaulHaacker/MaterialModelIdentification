function diff = diffFcn(par_norm,omega,storage_data,loss_data,weight_loss)
% difference of material model in frequency domain and data. for
% identification pruposes, namely to be handed to lsqnonlin

% Inputs:
%   par_norm      - Normalized parameters of the material model
%   omega         - frequencies at which the material properties are evaluated
%   storage_data  - Experimental storage modulus data in frequency domain
%   loss_data     - Experimental loss modulus data in frequency domain
%   weight_loss   - Weighting factor for the loss modulus data
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