% Test for G1StressDriven_DoubleOrderModelNonlin

clear
close all

alpha_1 = .3;
alpha_2 = .2;
E_0 = 100;
E_1 = 20;
E_2 = 20;
p_1 = 40;
p_2 = 40;
G = 20;
par =   [ alpha_1;
        alpha_2 ;
        E_0 ;
        E_1 ;
        E_2 ;
        p_1 ;
        p_2 ;
        G]; % parameters

time = 0:.0001:1; % time array
% stress = min(time / 0.01, 1); % ramp up to 1 at 0.5 seconds and stay constant
stress = ones(size(time)); 
% stress = zeros(size(time)); 
% stress(1:10)=1;
strain_0 = 0;

strain = G1StressDriven_DoubleOrderModelNonlin(par,stress,time,strain_0);

% plot results

figure
% tiledlayout('flow')
% nexttile
% plot(time,stress)
% xlabel('time $t$')
% ylabel('stress $\sigma(t)$')
% nexttile
plot(time,strain,'LineWidth',1.5)
xlabel('time $t$')
ylabel('strain $\varepsilon(t)$')
title('G1-algorithm applied to DOM')