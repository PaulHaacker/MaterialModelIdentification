function ComplexModulus = ComplexModSimulation_SingleOrderModelNonlin(par,omega)
% script finding Storage and Loss Modulus of the nl-SOM at given
% frequencies through simulation and discrete fourier transformation
% par  ...          (5-by-1)-array of parameters, where 
%                   alpha = par(1) \in (0,1)
%                   E0 = par(2)
%                   E1 = par(3)
%                   p1 = par(4)
%                   G = par(5)
% omega ...         (1xN) array of frequencies of exitation > 0

% extract length
N = length(omega);
ComplexModulus = zeros(size(omega));
for kk = 1:N
    dt = pi/omega(kk)/10;
    time = 0:dt:10^4*dt;  % Time vector 
    strain = @(t)  sin(omega(kk)* t);  % sinusoidal strain function
    stress_0 = 0;  % Initial stress condition
    % find numerical solution of nonlinear SOM
    [time_out, stress_vec] = G1StrainDriven_SingleOrderModelNonlin(par, strain, time, stress_0);

    T = 2*pi/omega(kk);
    x_samples = stress_vec(end-T/dt:end);
    % complex modulus are the first fourier coefficents
    ComplexModulus(kk)= (2/T) * dt*trapz(x_samples .*(-1i).* exp(1i*omega(kk) * time_out(end-T/dt:end)));
end
end