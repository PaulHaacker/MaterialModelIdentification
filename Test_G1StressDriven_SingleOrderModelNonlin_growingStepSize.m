clear
close all

alpha = 1/3;
E_0 = 2;
E_1 = 1;
p_1 = 0.5;
G=1;
par =   [alpha;
        E_0;
        E_1;
        p_1;
        G]; % parameters

tspan = [0,1]; % time span
stress = @(t) t.^0;
strain_0 = 0;

[t_vec,strain_vec] = G1StressDriven_SingleOrderModelNonlin_growingStepSize(par,stress,tspan,strain_0);

% plot results

figure
plot(t_vec,strain_vec,'o-')
xlabel('time $t$')
ylabel('strain $\varepsilon(t)$')
title('G1-algorithm applied to $D^\alpha \sigma + b\sigma = cD^\alpha \varepsilon + d\varepsilon$')