function [c,ceq] = nonlincon_SingleOrderModel(p)
% nonlinear constraint to the normalized parameters b, c, d for them to
% correspond to physical parameters.
c = p(4)-p(3)*p(2);     % Compute nonlinear inequalities at p.
ceq = [];   % Compute nonlinear equalities at p.
end