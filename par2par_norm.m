function par_norm = par2par_norm(par)
% transforming physical parameters to their normalized representation of
% the complex modulus
alpha = par(1);
E_0 = par(2);
E_1 = par(3);
p_1 = par(4);
par_norm =   [alpha;
            E_1/p_1;
            E_0 + E_1;
            E_0*E_1/p_1];
end