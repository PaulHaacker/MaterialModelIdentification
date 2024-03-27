function par_norm_identified = identify_MaterialModel(omega_data, storage_data, loss_data)
    % Define parameter bounds
    lb = zeros(4,1);
    ub = [1; Inf; Inf; Inf];

    % No linear equality or inequality constraints
    A_ineq = [];
    b_ineq = [];
    Aeq = [];
    beq = [];

    % Initial guess
    par_0 = ones(4,1);
    par_norm0 = par2par_norm(par_0); % Convert to normalized parameters

    % Weight for loss modulus
    weight_loss = 1;

    % Options for lsqnonlin
    options = optimoptions("lsqnonlin","Algorithm","trust-region-reflective");

    % Call lsqnonlin to identify parameters
    par_norm_identified = lsqnonlin(@(p) diffFcn(p, omega_data, storage_data, loss_data, weight_loss), ...
        par_norm0, lb, ub, A_ineq, b_ineq, Aeq, beq, @(p) nonlincon_identification(p), options);
end