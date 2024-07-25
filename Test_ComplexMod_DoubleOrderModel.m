%% Test call of ComplexModulusFcn
clear
close all

% % % Note: parameters E_0, E_1, p_1 are unconstrained on > 0, but this
% implies that parameters b, c, d must satisfy c-d/b > 0 ! Else the model
% is unphysical!
alpha_1 = .5;
alpha_2 = 1/2;
E_0 = 100;
E_1 = 20;
E_2 = 20;
p_1 = 40;
p_2 = 40;
par =   [ alpha_1;
        alpha_2 ;
        E_0 ;
        E_1 ;
        E_2 ;
        p_1 ;
        p_2 ]; % parameters

omega = logspace(-10,10,100); % frequency range
ComplexModulus = ComplexMod_SingleOrderModel(par,omega);

% plot
figure
tiledlayout('flow')
nexttile
plot(omega,real(ComplexModulus),'o-',omega,imag(ComplexModulus),'o-');
set(gca,'xscale','log')
grid on
legend('storage modulus $E''= \Re\{E^\ast\}$','loss modulus $E''''= \Im\{E^\ast\}$'  )
xlabel('frequency in ???')