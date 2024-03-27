%% Identification of fractional material model from DMA (dynamic mechanical analysis),
% experiments conducted by Ondrej Farkas and Alexander Lion (Uni-BW Munich)
clear
close all

% % load data
% load('masterData20deg.mat')

% generate "data" directly from the model to test the identification - this
% is essentially just a copy of the test script of the model.
alpha = 1/3;
E_0 = 2;
E_1 = 1;
p_1 = 0.5;
par_test =   [alpha;
            E_0;
            E_1;
            p_1]; % parameters in original form
par_norm_test = par2par_norm(par_test);

omega_data = logspace(-5,5,100); % frequency range
ComplexModulus_data = ComplexModulusFcn(par_norm_test,omega_data);
storage_data = real(ComplexModulus_data);
loss_data = imag(ComplexModulus_data);

% parameters to be identified are the normalized parameters which are an
% input to the model "ComplexModulusFcn" and live in R^4, namely
% par_norm ...      (4-by-1)-array of normalized parameters, where 
%                   alpha = par_norm(1) \in (0,1)
%                   b = par_norm(2) = E_1/p_1 > 0
%                   c = par_norm(3) = E_0 + E_1 > 0
%                   d = par_norm(4) = E_0*E_1/p_1> 0
% The parameters are subject to constraints, namely they are partially
% bounded:
lb = zeros(4,1);
ub = [1 inf inf inf]';

% no linear equality or inequality constraints
A_ineq = [];
b_ineq = [];
Aeq = [];
beq = [];

% nonlinear inequality constraint is defined in "nonlincon_identification"

% initial guess
par_0 = ones(4,1);
par_norm0 = par2par_norm(par_0); % normalized parameters

% weight for loss modulus - needed since magnitude of loss data generally
% is much smaller than storage data.
weight_loss = 1;

options = optimoptions("lsqnonlin","Algorithm","trust-region-reflective");
par_norm_lsqnonlin = lsqnonlin(@(p)diffFcn(p,omega_data,storage_data,loss_data,weight_loss) ...
    ,par_norm0,lb,ub,A_ineq,b_ineq,Aeq,beq,...
    @(p) nonlincon_identification(p),options);

