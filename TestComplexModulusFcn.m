%% Test call of ComplexModulusFcn
clear
close all

% documentation of the function to be tested:
% Material Model of a Viscoelastic Material in the frequency domain
% inputs:
% par_norm ...      (4-by-1)-array of normalized parameters, where 
%                   alpha = par_norm(1) \in (0,1)
%                   b = par_norm(2) = E_1/p_1 > 0
%                   c = par_norm(3) = E_0 + E_1 > 0
%                   d = par_norm(4) = E_0*E_1/p_1> 0
%
% omega ...         frequency of exitation > 0

% % % Note: parameters E_0, E_1, p_1 are unconstrained on > 0, but this
% implies that parameters b, c, d must satisfy c-d/b > 0 ! Else the model
% is unphysical!
alpha = 1/3;
E_0 = 2;
E_1 = 1;
p_1 = 0.5;
par =   [alpha;
        E_0;
        E_1;
        p_1]; % parameters

omega = logspace(-5,5,100); % frequency range
ComplexModulus = ComplexModulusFcn(par2par_norm(par),omega);

% plot
figure
tiledlayout('flow')
nexttile
plot(omega,real(ComplexModulus),'o-',omega,imag(ComplexModulus),'o-');
set(gca,'xscale','log')
grid on
legend('storage modulus $E"= \Re\{E^\ast\}$','loss modulus $E""= \Im\{E^\ast\}$'  )
xlabel('frequency in ???')