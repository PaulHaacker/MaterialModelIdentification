function stress = G1StrainDriven_DoubleOrderModelNonlin(par,strain,time,stress_0)
% Double Order Material Model of a Viscoelastic Material, time-domain
% response using G1-algorithm - Strain Driven. Solving for Stress Answer.
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
% strain ...        array of strain signal, 1-by-N_total
% stress_0 ...      initial condition of stress

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
if length(strain)~= N_total
    error('invalid input length')
end
dt = diff(time);
tolerance = 1e-10; % Define a small tolerance value
is_equidistant = all(abs(dt - dt(1)) < tolerance);
if ~is_equidistant
    error('invalid time array')
end
dt = dt(1); % time steps

stress = zeros(size(strain));
stress(1) = stress_0;

A_Grunwald_1_mod = zeros(size(stress));
A_Grunwald_1_mod(1) = dt^(-alpha_1);
A_Grunwald_2_mod = zeros(size(stress));
A_Grunwald_2_mod(1) = dt^(-alpha_2);
A_Grunwald_3_mod = zeros(size(stress));
A_Grunwald_3_mod(1) = dt^(-alpha_3);

% time stepping
for kk = 2:N_total
    % compute Grunwald coefficients recursively
    A_Grunwald_1_mod(kk)= A_Grunwald_1_mod(kk-1)*(kk-2-alpha_1)/(kk-1);
    A_Grunwald_2_mod(kk)= A_Grunwald_2_mod(kk-1)*(kk-2-alpha_2)/(kk-1);
    A_Grunwald_3_mod(kk)= A_Grunwald_3_mod(kk-1)*(kk-2-alpha_3)/(kk-1);
    % compute fractional derivative terms of stress
    der_stress = sum((a10*A_Grunwald_1_mod(2:kk)...
        +a01*A_Grunwald_2_mod(2:kk)...
        +a11*A_Grunwald_3_mod(2:kk)).*fliplr(stress(1:kk-1)));
    % compute a part of the fractional der. of strain
    der_term_strain = sum( (b10*A_Grunwald_1_mod(1:kk)...
        +b01*A_Grunwald_2_mod(1:kk)...
        +b11*A_Grunwald_3_mod(1:kk)).*fliplr(strain(1:kk)) ...
        +G*(a10*A_Grunwald_1_mod(1:kk)...
        +a01*A_Grunwald_2_mod(1:kk)...
        +a11*A_Grunwald_3_mod(1:kk)).*fliplr(strain(1:kk).^3) );
    
    % der_term_strain = dt^(-alpha)*sum(A_Grunwald_1(2:kk).*fliplr(strain(1:kk-1)));
    % time step: find next value of strain
    stress(kk) = (-der_stress+a0*(E_0*strain(kk)+G*strain(kk)^3)+der_term_strain)...
        /(a0+a11*dt^-alpha_3+a01*dt^-alpha_2+a10*dt^-alpha_1);
end

end