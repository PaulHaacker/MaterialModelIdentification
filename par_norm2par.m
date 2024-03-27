function par = par_norm2par(par_norm)
% transforming physical parameters to their normalized representation of
% the complex modulus
alpha = par_norm(1);
b = par_norm(2);
c = par_norm(3);
d = par_norm(4);
par =   [alpha;
            d/b;
            c-d/b;
            (c-d/b)/b]; 
end