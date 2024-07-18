function [t_vec,strain_vec] = G1StressDriven_SingleOrderModel_growingStepSize(par_norm,stress,tspan,strain_0)
% Single Order Material Model of a Viscoelastic Material, time-domain
% response using G1-algorithm
% inputs:
% par_norm ...      (4-by-1)-array of normalized parameters, where 
%                   alpha = par_norm(1) \in (0,1)
%                   b = par_norm(2) = E_1/p_1 > 0
%                   c = par_norm(3) = E_0 + E_1 > 0
%                   d = par_norm(4) = E_0*E_1/p_1> 0
% tspan ...         time span of the form [t_0 t_1]
% stress ...        function of stress signal, takes inputs from t_0 to t_1
% strain_0 ...      initial condition of strain (scalar)

% % % Note: parameters E_0, E_1, p_1 are unconstrained on > 0, but this
% implies that parameters b, c, d must satisfy bc-d > 0 ! Else the model
% is not physically sound!

% extract parameters
alpha = par_norm(1);
b = par_norm(2);
c = par_norm(3);
d = par_norm(4);

% paramters of time stepping
N_c = 20; % number of steps performed with constant dt
N_d = 10; % number of times the time step is doubled and the constant time step integration is performed, = N_d - 1

T = tspan(2)-tspan(1);

dt = T/N_c/sum((2*ones(1,N_d)).^(0:N_d-1));
t_vec = zeros(N_c*N_d+1,1);
t_vec(1) = tspan(1);

strain_vec = zeros(size(t_vec));
strain_vec(1) = strain_0;
strain_inner = strain_0;
strain_hist = strain_0;

% Generate t_vec
index = 2;
current_time = tspan(1);
current_dt = dt;

for ii = 1:N_d
    new_times = current_time + (1:N_c) * current_dt;
    t_vec(index:index+N_c-1) = new_times;
    index = index + N_c;
    current_time = new_times(end);
    current_dt = current_dt * 2;
end

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

% time stepping
for kk = 1:N_d
    stress_inner = [stress_hist;stress_vec(2+N_c*(kk-1):1+N_c*kk)];
    strain_inner = [strain_hist;zeros(N_c,1)];

    for jj = 1:N_c
        % compute fractional derivative of stress
        length_hist_curr = length(stress_hist)+jj;
        der_stress = dt^(-alpha)*sum(A_Grunwald(1:length_hist_curr).*flip(stress_inner(1:length_hist_curr)));
        % compute a part of the fractional der. of strain
        der_term_strain = dt^(-alpha)*sum(A_Grunwald(2:length_hist_curr).*flip(strain_inner(1:length_hist_curr-1)));
        % time step: find next value of strain
        strain_inner(length_hist_curr) = (der_stress+b*stress_inner(length_hist_curr)-c*der_term_strain)/(d+c*dt^(-alpha));
    end
    stress_hist = every2ndentry(stress_inner);
    strain_hist = every2ndentry(strain_inner);
    strain_vec(2+N_c*(kk-1):1+N_c*kk) = strain_inner(end-N_c+1:end);
    dt = 2*dt;
end

end

function new = every2ndentry(old)
    help = fliplr(old);
    new = fliplr(help(1:2:end));
end