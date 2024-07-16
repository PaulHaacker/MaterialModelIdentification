function strain = G1StressDriven_SingleOrderModel(par_norm,stress,time,strain_0)
% Single Order Material Model of a Viscoelastic Material, time-domain
% response using G1-algorithm
% inputs:
% par_norm ...      (4-by-1)-array of normalized parameters, where 
%                   alpha = par_norm(1) \in (0,1)
%                   b = par_norm(2) = E_1/p_1 > 0
%                   c = par_norm(3) = E_0 + E_1 > 0
%                   d = par_norm(4) = E_0*E_1/p_1> 0
% time ...          array of time steps, 1-by-N_total. Must be equidistant.
% stress ...        array of stress signal, 1-by-N_total
% strain_0 ...      initial condition of strain

% % % Note: parameters E_0, E_1, p_1 are unconstrained on > 0, but this
% implies that parameters b, c, d must satisfy bc-d > 0 ! Else the model
% is not physically sound!

% extract parameters
alpha = par_norm(1);
b = par_norm(2);
c = par_norm(3);
d = par_norm(4);

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

A_Grunwald = zeros(size(strain));
A_Grunwald(1) = 1;

% time stepping
for kk = 2:N_total
    % compute Grunwald coefficients recursively
    A_Grunwald(kk)= A_Grunwald(kk-1)*(kk-2-alpha)/(kk-1);
    % compute fractional derivative of stress
    der_stress = dt^(-alpha)*sum(A_Grunwald(1:kk).*fliplr(stress(1:kk)));
    % compute a part of the fractional der. of strain
    der_term_strain = dt^(-alpha)*sum(A_Grunwald(2:kk).*fliplr(strain(1:kk-1)));
    % time step: find next value of strain
    strain(kk) = (der_stress+b*stress(kk)-c*der_term_strain)/(d+c*dt^(-alpha));
end

end