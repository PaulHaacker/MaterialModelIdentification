function strain = G1StressDriven_DoubleOrderModelNonlin(par,stress,time,strain_0)
% Single Order Material Model of a Viscoelastic Material, time-domain
% response using G1-algorithm
% inputs:
% par ...      (7-by-1)-array of normalized parameters, where 
%                   alpha_1 = par(1) \in (0,1)
%                   alpha_2 = par(2) \in (0,1)
%                   E_0 = par(3)
%                   E_1 = par(4)
%                   E_2 = par(5)
%                   p_1 = par(6)
%                   p_2 = par(7)
%                   G = par(8)
% time ...          array of time steps, 1-by-N_total. Must be equidistant.
% stress ...        array of stress signal, 1-by-N_total
% strain_0 ...      initial condition of strain

% extract parameters
alpha_1 = par(1);
alpha_2 = par(2);
E_0 = par(3);
E_1 = par(4);
E_2 = par(5);
p_1 = par(6);
p_2 = par(7);
G = par(8);

alpha_3 = alpha_1+alpha_2;

a11 = 1/E_1/E_2;
a01 = 1/E_2/p_1;
a10 = 1/E_1/p_2;
a0 = 1/p_1/p_2;
b11= 1/E_1 + 1/E_2 + E_0/E_1/E_2;
b01 = (1+E_0/E_2)/p_1;
b10 = (1+E_0/E_1)/p_2;
b0 = E_0/p_1/p_2;

N_total = length(time);
if length(stress)~= N_total
    error('invalid input length')
end
dt = diff(time);
tolerance = 1e-10; % Define a small tolerance value
is_equidistant = all(abs(dt - dt(1)) < tolerance);
if ~is_equidistant
    error('invalid time array')
end
dt = dt(1); % time steps

strain = zeros(size(stress));
strain(1) = strain_0;

A_Grunwald_1_tilde = zeros(size(strain));
A_Grunwald_1_tilde(1) = dt^-alpha_1;
A_Grunwald_2_tilde = zeros(size(strain));
A_Grunwald_2_tilde(1) = dt^-alpha_2;
A_Grunwald_3_tilde = zeros(size(strain));
A_Grunwald_3_tilde(1) = dt^-alpha_3;

% time stepping
for kk = 2:N_total
    % compute Grunwald coefficients recursively
    A_Grunwald_1_tilde(kk)= A_Grunwald_1_tilde(kk-1)*(kk-2-alpha_1)/(kk-1);
    A_Grunwald_2_tilde(kk)= A_Grunwald_2_tilde(kk-1)*(kk-2-alpha_2)/(kk-1);
    A_Grunwald_3_tilde(kk)= A_Grunwald_3_tilde(kk-1)*(kk-2-alpha_3)/(kk-1);
    % compute fractional derivative terms of stress
    der_stress = sum((a10*A_Grunwald_1_tilde(2:kk)...
        +a01*A_Grunwald_2_tilde(2:kk)...
        +a11*A_Grunwald_3_tilde(2:kk)).*fliplr(stress(1:kk-1)));
    instant_stress = (a0+a11*dt^-alpha_3+a01*dt^-alpha_2+a10*dt^-alpha_1)*stress(kk);
    % compute a part of the fractional der. of strain
    der_term_strain = sum((b10*A_Grunwald_1_tilde(2:kk)...
        +b01*A_Grunwald_2_tilde(2:kk)...
        +b11*A_Grunwald_3_tilde(2:kk)).*fliplr(strain(1:kk-1)) ...
        +G*(a10*A_Grunwald_1_tilde(2:kk)...
        +a01*A_Grunwald_2_tilde(2:kk)...
        +a11*A_Grunwald_3_tilde(2:kk)).*fliplr(strain(1:kk-1).^3));
    
    % der_term_strain = dt^(-alpha)*sum(A_Grunwald_1(2:kk).*fliplr(strain(1:kk-1)));
    % time step: find next value of strain
    df = @(x)(b0+b11*dt^-alpha_3+b01*dt^-alpha_2+b10*dt^-alpha_1) ...
        + 2*(a0+a11*dt^-alpha_3+a01*dt^-alpha_2+a10*dt^-alpha_1)*x.^2 ;
    f = @(x) (b0+b11*dt^-alpha_3+b01*dt^-alpha_2+b10*dt^-alpha_1)*x ...
        + (a0+a11*dt^-alpha_3+a01*dt^-alpha_2+a10*dt^-alpha_1)*x.^3 ...
        + der_term_strain - der_stress - instant_stress;
    strain(kk) = newtonRaphson(f, df, strain(kk-1), 10^-6, 30);

    % strain(kk) = (der_stress+a0*stress(kk)-der_term_strain)...
    %     /(b0+b11*dt^-alpha_3+b01*dt^-alpha_2+b10*dt^-alpha_1);
end

end


function root = newtonRaphson(f, df, x0, tol, maxIter)
    % f - The function for which we are finding the root
    % df - The derivative of the function f
    % x0 - Initial guess
    % tol - Tolerance for convergence
    % maxIter - Maximum number of iterations allowed
    
    % Initialize the iteration counter
    iter = 0;
    % Set the initial guess as the current guess
    x = x0;
    
    % Iterate until convergence or max iterations
    while iter < maxIter
        % Evaluate the function and its derivative at the current guess
        fx = f(x);
        dfx = df(x);
        
        % Newton-Raphson update step
        x_new = x - fx / dfx;
        
        % Check for convergence
        if abs(x_new - x) < tol
            root = x_new;
            % fprintf('Converged to root: %.6f\n', root);
            return;
        end
        
        % Update the guess for the next iteration
        x = x_new;
        
        % Increment the iteration counter
        iter = iter + 1;
    end
    
    % If we reach max iterations without convergence
    warning('Maximum iterations reached without convergence.');
    root = x;
end
