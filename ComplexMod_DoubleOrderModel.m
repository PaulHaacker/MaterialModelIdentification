function ComplexModulus = ComplexMod_DoubleOrderModel(par,omega)
% Material Model of a Viscoelastic Material in the frequency domain
% inputs:
% par ...      (7-by-1)-array of normalized parameters, where 
%                   alpha_1 = par(1) \in (0,1)
%                   alpha_2 = par(2) \in (0,1)
%                   E_0 = par(3)
%                   E_1 = par(4)
%                   E_2 = par(5)
%                   p_1 = par(6)
%                   p_2 = par(7)
%
% omega ...         frequency of exitation > 0

alpha_1 = par(1);
alpha_2 = par(2);
E_0 = par(3);
E_1 = par(4);
E_2 = par(5);
p_1 = par(6);
p_2 = par(7);

a11 = 1/E_1/E_2;
a01 = 1/E_1/p_2;
a10 = 1/E_2/p_1;
a0 = 1/p_1/p_2;
b11= 1/E_1 + 1/E_2 + E_0/E_1/E_2;
b01 = (1+E_0/E_2)/p_1;
b10 = (1+E_0/E_1)/p_2;
b0 = E_0/p_1/p_2;

s = 1j*omega;

ComplexModulus = (b11*s.^(alpha_1+alpha_2)+b01*s.^alpha_2+b10*s.^alpha_1+b0)./...
    (a11*s.^(alpha_1+alpha_2)+a01*s.^alpha_2+a10*s.^alpha_1+a0);
end