function [time,stress_vec] = G1StrainDriven_SingleOrderModelNonlin(par,strain,time,stress_0)
% Single Order nonlinear Material Model (nl-SOM) of a Viscoelastic Material, time-domain
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
% t_vec ... 	    (1-by-N)-array of time instances
% stress_vec ...    (1-by-N)-array of simulated stress values at corresponding time instances

%% extract parameters + initializing
alpha = par(1);
E0 = par(2);
E1 = par(3);
p1 = par(4);
G = par(5);

% paramters of time stepping
N_total = length(time); % number of steps performed with constant dt

h = diff(time);
tolerance = 1e-6; % Define a small tolerance value
is_equidistant = all(abs(h - h(1)) < tolerance);
if ~is_equidistant
    error('invalid time array')
end
h = h(1); % time steps

stress_vec = zeros(size(time));
stress_vec(1) = stress_0;

strain_vec = strain(time);

%% run G1

A_Grunwald = zeros(size(time));
A_Grunwald(1) = 1;

help_b = h^alpha*E0*E1/p1 + E0 +E1;
help_c = h^alpha*E1/p1+1;

% time stepping
for kk = 2:N_total
    % compute Grunwald coefficients recursively
    A_Grunwald(kk)= A_Grunwald(kk-1)*(kk-2-alpha)/(kk-1);
    % compute fractional derivative terms
    der_terms = sum(A_Grunwald(2:kk).*fliplr((E0+E1)*strain_vec(2:kk)+G*strain_vec(2:kk).^3-stress_vec(2:kk)));
    % time step: find next value of strain
    stress_vec(kk) = (der_terms+help_b*strain_vec(kk)+help_c*strain_vec(kk)^3)/help_c;
end

end

% function new = every2ndentry(old)
%     help = fliplr(old);
%     new = fliplr(help(1:2:end));
% end
% 
% function root = newtonRaphson(f, df, x0, tol, maxIter)
%     % f - The function for which we are finding the root
%     % df - The derivative of the function f
%     % x0 - Initial guess
%     % tol - Tolerance for convergence
%     % maxIter - Maximum number of iterations allowed
% 
%     % Initialize the iteration counter
%     iter = 0;
%     % Set the initial guess as the current guess
%     x = x0;
% 
%     % Iterate until convergence or max iterations
%     while iter < maxIter
%         % Evaluate the function and its derivative at the current guess
%         fx = f(x);
%         dfx = df(x);
% 
%         % Newton-Raphson update step
%         x_new = x - fx / dfx;
% 
%         % Check for convergence
%         if abs(x_new - x) < tol
%             root = x_new;
%             % fprintf('Converged to root: %.6f\n', root);
%             return;
%         end
% 
%         % Update the guess for the next iteration
%         x = x_new;
% 
%         % Increment the iteration counter
%         iter = iter + 1;
%     end
% 
%     % If we reach max iterations without convergence
%     warning('Maximum iterations reached without convergence.');
%     root = x;
% end
