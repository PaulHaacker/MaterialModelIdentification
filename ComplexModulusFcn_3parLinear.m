function ComplexModulus = ComplexModulusFcn_3parLinear(par_norm,omega)
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

% extract parameters
alpha = par_norm(1);
b = par_norm(2);
c = par_norm(3);
d = par_norm(4);

ComplexModulus = (c*(1j*omega).^alpha + d)./((1j*omega).^alpha + b);
end