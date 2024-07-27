function strain = G1StressDriven_DoubleOrderModel(par,stress,time,strain_0)
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

alpha_3 = alpha_1+alpha_2;

a11 = 1/E_1/E_2;
a01 = 1/E_1/p_2;
a10 = 1/E_2/p_1;
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

A_Grunwald_1 = zeros(size(strain));
A_Grunwald_1(1) = 1;
A_Grunwald_2 = zeros(size(strain));
A_Grunwald_2(1) = 1;
A_Grunwald_3 = zeros(size(strain));
A_Grunwald_3(1) = 1;

% time stepping
for kk = 2:N_total
    % compute Grunwald coefficients recursively
    A_Grunwald_1(kk)= A_Grunwald_1(kk-1)*(kk-2-alpha_1)/(kk-1);
    A_Grunwald_2(kk)= A_Grunwald_2(kk-1)*(kk-2-alpha_2)/(kk-1);
    A_Grunwald_3(kk)= A_Grunwald_3(kk-1)*(kk-2-alpha_3)/(kk-1);
    % compute fractional derivative terms of stress
    der_stress = sum((a10*dt^(-alpha_1)*A_Grunwald_1(1:kk)...
        +a01*dt^(-alpha_2)*A_Grunwald_2(1:kk)...
        +a11*dt^(-alpha_3)*A_Grunwald_3(1:kk)).*fliplr(stress(1:kk)));
    % compute a part of the fractional der. of strain
    der_term_strain = sum((b10*dt^(-alpha_1)*A_Grunwald_1(2:kk)...
        +b01*dt^(-alpha_2)*A_Grunwald_2(2:kk)...
        +b11*dt^(-alpha_3)*A_Grunwald_3(2:kk)).*fliplr(strain(1:kk-1)));
    
    % der_term_strain = dt^(-alpha)*sum(A_Grunwald_1(2:kk).*fliplr(strain(1:kk-1)));
    % time step: find next value of strain
    strain(kk) = (der_stress+a0*stress(kk)-der_term_strain)...
        /(b0+b11*dt^-alpha_3+b01*dt^-alpha_2+b10*dt^-alpha_1);
end

end