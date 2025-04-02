clc;
clear;
close all;

% Define test parameters
alpha = .6;
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
semilogx(omega, abs(a),'o-');
hold on
plot(omega,abs(real(ComplexModulus)),'o-');
xlabel('Frequency $\omega$ (rad/s)', 'FontSize', 14);
ylabel('Storage Modulus', 'FontSize', 14);
legend('harmonic balance $|a(\omega)|$', 'Fourier transform $\mathrm{Re}\bar E$', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 14);

subplot(2,1,2);
semilogx(omega, abs(b), 'o-');
hold on
plot(omega,abs(imag(ComplexModulus)),'o-');
xlabel('Frequency $\omega$ (rad/s)', 'FontSize', 14);
ylabel('Loss Modulus', 'FontSize', 14);
legend('harmonic balance $|b(\omega)|$', 'Fourier transform $\mathrm{Im}\bar E$', 'FontSize', 14);
grid on;

% Set font size for all axes
set(gca, 'FontSize', 14);


% disp('Test completed successfully.');
