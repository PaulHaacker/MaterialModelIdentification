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

time = 0:.001:1; % time array
stress = ones(size(time));
strain_0 = 0;

strain = G1StressDriven_SingleOrderModel(par2par_norm(par),stress,time,strain_0);

% plot results

figure
plot(time,strain)
xlabel('time $t$')
ylabel('strain $\varepsilon(t)$')
title('G1-algorithm applied to $D^\alpha \sigma + b\sigma = cD^\alpha \varepsilon + d\varepsilon$')