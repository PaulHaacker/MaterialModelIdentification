function [t_vec,stress_vec] = G1StrainDriven_SingleOrderModelNonlin(par,strain,time,stress_0)
% Single Order Material Model of a Viscoelastic Material, time-domain
% response using G1-algorithm
% inputs:
% par_     ...      (5-by-1)-array of parameters, where 
%                   alpha = par(1) \in (0,1)
%                   E0 = par(2)
%                   E1 = par(3)
%                   p1 = par(4)
%                   G = par(5)
% time ...          array of time steps, 1-by-N_total. Must be equidistant.
% strain ...        function of strain signal, takes inputs from t_0 to t_1
% stress_0 ...      initial condition of stress (scalar)

% outputs:
% t_vec ... 	    (1-by-N)-array of time instances, where the number N is hard-coded within the function
% stress_vec ...    (1-by-N)-array of simulated stress values at corresponding time instances

% extract parameters
alpha = par(1);
E0 = par(2);
E1 = par(3);
p1 = par(4);
G = par(5);

% paramters of time stepping
N_c = length(time); % number of steps performed with constant dt
N_d = 1; % number of times the time step is doubled and the constant time step integration is performed, = N_d - 1

dt = diff(time);
tolerance = 1e-10; % Define a small tolerance value
is_equidistant = all(abs(dt - dt(1)) < tolerance);
if ~is_equidistant
    error('invalid time array')
end
dt = dt(1); % time steps

t_vec = time;

strain_vec = zeros(size(t_vec));
strain_vec(1) = strain_0;
strain_inner = strain_0;
strain_hist = strain_0;

% % Generate t_vec
% index = 2;
% current_time = time(1);
% current_dt = dt;
% 
% for ii = 1:N_d
%     new_times = current_time + (1:N_c) * current_dt;
%     t_vec(index:index+N_c-1) = new_times;
%     index = index + N_c;
%     current_time = new_times(end);
%     current_dt = current_dt * 2;
% end

stress_vec = stress(t_vec);
stress_hist = stress_vec(1);

% % find number of Grunwald coeffs necessary
% sum_G = 2;
% for kk =1:N_d
%     sum_G = floor(0.5*sum_G) + N_c;
% end
sum_G = N_c*N_d+1;

% compute Grunwald coeffs
A_Grunwald = zeros(sum_G,1);
A_Grunwald(1) = 1;

for kk = 2:sum_G
    % compute Grunwald coefficients recursively
    A_Grunwald(kk)= A_Grunwald(kk-1)*(kk-2-alpha)/(kk-1);
end

options = optimset('Display','off');

% time stepping
for kk = 1:N_d
    stress_inner = stress_vec';
    strain_inner = [strain_hist;zeros(N_c-1,1)];

    for jj = 1:N_c-1
        % compute fractional derivative of stress
        length_hist_curr = length(stress_hist)+jj;
        der_stress = dt^(-alpha)*sum(A_Grunwald(1:length_hist_curr).*flip(stress_inner(1:length_hist_curr)));
        % compute a part of the fractional der. of strain
        der_term_strain = sum(A_Grunwald(2:length_hist_curr).*...
            flip((E0+E1)*strain_inner(1:length_hist_curr-1) + G*strain_inner(1:length_hist_curr-1).^3));
        
        % % time step: find next value of strain

        % strain_inner(length_hist_curr) = fsolve(...
        %     @(x)der_term_strain - dt^alpha*(der_stress + E1/p1*stress_inner(length_hist_curr))...
        %     + (E0+E1+dt^alpha*E0*E1/p1)*x+ G*(1+dt^alpha*E1/p1)*x^3, ...
        %     strain_inner(length_hist_curr-1),options);

        % strain_inner(length_hist_curr) = fzero(...
        %     @(x)der_term_strain - dt^alpha*(der_stress + E1/p1*stress_inner(length_hist_curr))...
        %     + (E0+E1+dt^alpha*E0*E1/p1)*x+ G*(1+dt^alpha*E1/p1)*x.^3, ...
        %     strain_inner(length_hist_curr-1),options);
        
        df = @(x)(E0+E1+dt^alpha*E0*E1/p1)+ 3*G*(1+dt^alpha*E1/p1)*x.^2;
        f = @(x)der_term_strain - dt^alpha*(der_stress + E1/p1*stress_inner(length_hist_curr))...
            + (E0+E1+dt^alpha*E0*E1/p1)*x+ G*(1+dt^alpha*E1/p1)*x.^3;
        strain_inner(length_hist_curr) = newtonRaphson(f, df, strain_inner(length_hist_curr-1), 10^-6, 30);

        % p = [G*(1+dt^alpha*E1/p1),0,(E0+E1+dt^alpha*E0*E1/p1)...
        %     ,der_term_strain - dt^alpha*(der_stress + E1/p1*stress_inner(length_hist_curr))];
        % all_roots = roots(p);
        % [~, idx] = min(abs(all_roots - strain_inner(length_hist_curr-1))); % Find index of the closest root
        % strain_inner(length_hist_curr) = all_roots(idx); % Select the corresponding root
    end
    stress_hist = every2ndentry(stress_inner);
    strain_hist = every2ndentry(strain_inner);
    strain_vec = strain_inner;
    dt = 2*dt;
end

end

function new = every2ndentry(old)
    help = fliplr(old);
    new = fliplr(help(1:2:end));
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
