clc;
clear;
close all;

% Define test parameters
alpha = .4;
E_0 = 3;
E_1 = 1;
p_1 = 1;
par_NL =   [alpha;
        E_0;
        E_1;
        p_1;
        0]; % parameters
par_lin =   [alpha;
        E_0;
        E_1;
        p_1]; % parameters

omega = logspace(-5,10,100); % frequency range

% Call the function
[a, b] = HB_nlSOM(par_NL, omega);
ComplexModulus = ComplexMod_SingleOrderModel(par2par_norm(par_lin),omega);

% % Display results
% disp('Sine Coefficients (Storage Modulus):');
% disp(a);
% 
% disp('Cosine Coefficients (Loss Modulus):');
% disp(b);

% % Validate outputs
% assert(isvector(a) && length(a) == length(omega), 'Output "a" should be a vector of the same size as omega.');
% assert(isvector(b) && length(b) == length(omega), 'Output "b" should be a vector of the same size as omega.');
% assert(all(isfinite(a)), 'Output "a" contains non-finite values.');
% assert(all(isfinite(b)), 'Output "b" contains non-finite values.');

% Plot the results
figure;
subplot(2,1,1);
semilogx(omega, a, 'bo-', 'LineWidth', 1.5);
hold on
plot(omega,real(ComplexModulus),'o-');
xlabel('Frequency $\omega$ (rad/s)');
ylabel('Storage Modulus $a(\omega)$');
legend('harmonic balance','transfer function')
title('Storage Modulus');
grid on;

subplot(2,1,2);
semilogx(omega, b, 'ro-', 'LineWidth', 1.5);
hold on
plot(omega,imag(ComplexModulus),'o-');
xlabel('Frequency $\omega$ (rad/s)');
ylabel('Loss Modulus $b(\omega)$');
legend('harmonic balance','transfer function')
title('Loss Modulus');
grid on;

% disp('Test completed successfully.');
