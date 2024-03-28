function [c,ceq] = nonlincon_3parLinear_alpha1(p)
% nonlinear constraint to the normalized parameters b, c, d for them to
% correspond to physical parameters.
c = p(4)/p(2)-p(3);     % Compute nonlinear inequalities at p.
ceq = p(1)-1;   % fix alpha = 1 to enforce whole-order model
end