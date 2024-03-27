function [c,ceq] = nonlincon_3parLinear(p)
% nonlinear constraint to the normalized parameters b, c, d for them to
% correspond to physical parameters.
c = p(4)/p(2)-p(3);     % Compute nonlinear inequalities at p.
ceq = [];   % Compute nonlinear equalities at p.
end