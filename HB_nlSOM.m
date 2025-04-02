function [a,b] = HB_nlSOM(par,omega)
    % Harmonic Balance of the nlSOM, assuming input \sigma(t) = sin(\omega
    % t). 
    % inputs:
    % par     ...      (5-by-1)-array of parameters, where 
    %                   alpha = par(1) \in (0,1)
    %                   E0 = par(2)
    %                   E1 = par(3)
    %                   p1 = par(4)
    %                   G = par(5)
    % omega   ...      (N-by-1)-array of radial frequencies
    %
    % outputs:
    % a       ...       (N-by-1)-array of sine coeff (storage modulus)
    % b       ...       (N-by-1)-array of cos coeff (loss modulus)
    
    alpha = par(1);
    E0 = par(2);
    E1 = par(3);
    p1 = par(4);
    G = par(5);
    N=length(omega);
    a = zeros(size(omega));
    b = zeros(size(omega));
    for kk = 1:N
        A_mat = [E1/p1+omega(kk)^alpha*cos(alpha*pi/2), -omega(kk)^alpha*sin(alpha*pi/2)
                omega(kk)^alpha*sin(alpha*pi/2),E1/p1+omega(kk)^alpha*cos(alpha*pi/2)];
        b_mat = [(E0+E1+3/4*G)*omega(kk)^alpha*cos(alpha*pi/2)+E0*E1/p1+3/4*G;
                (E0+E1+3/4*G)*omega(kk)^alpha*sin(alpha*pi/2)];
        soln = A_mat\b_mat;
        a(kk) = soln(1);
        b(kk) = soln(2);
    end

end

% for stress-driven:

% function [a,b] = HB_nlSOM(par,omega)
%     % Harmonic Balance of the nlSOM, assuming input \sigma(t) = sin(\omega
%     % t). 
%     % inputs:
%     % par     ...      (5-by-1)-array of parameters, where 
%     %                   alpha = par(1) \in (0,1)
%     %                   E0 = par(2)
%     %                   E1 = par(3)
%     %                   p1 = par(4)
%     %                   G = par(5)
%     % omega   ...      (N-by-1)-array of radial frequencies
%     %
%     % outputs:
%     % a       ...       (N-by-1)-array of sine coeff (storage modulus)
%     % b       ...       (N-by-1)-array of cos coeff (loss modulus)
% 
%     N=length(omega);
%     a = zeros(size(omega));
%     b = zeros(size(omega));
%     a0 = 1;
%     b0 = 1;
%     for kk = 1:N
%         [soln,res] = fsolve(@(coeffs)residual(coeffs,par,omega(kk)),[a0;b0]);
%         if norm(res)>10^(-4)
%             error('fsolve failed to solve balance equation')
%         end
%         a(kk) = soln(1);
%         b(kk) = soln(2);
%         a0 = soln(1);
%         b0 = soln(2);
%     end
% 
% end

% function r = residual(coeffs,par,omega)
%     r = zeros(2,1);
% 
%     a = coeffs(1);
%     b = coeffs(2);
% 
%     % extract parameters
%     alpha = par(1);
%     E0 = par(2);
%     E1 = par(3);
%     p1 = par(4);
%     G = par(5);
% 
%     alpha_tild = alpha*pi/2;
% 
%     r(1)= a*(E0*E1/p1 + (E0+E1)*omega^alpha*cos(alpha_tild))...
%         - b*omega^alpha*(E0+E1)*sin(alpha_tild)...
%         +(a^3+a*b^2)*3/4*G*(E1/p1+omega^alpha*cos(alpha_tild))...
%         -(b^3+a^2*b)*3/4*G*omega^alpha*sin(alpha_tild)...
%         -(omega^alpha*cos(alpha_tild)+E1/p1);
%     r(2)= b*(E0*E1/p1 + (E0+E1)*omega^alpha*sin(alpha_tild))...
%         +a*omega^alpha*(E0+E1)*sin(alpha_tild)...
%         +(b^3+a^2*b)*3/4*G*(E1/p1+omega^alpha*cos(alpha_tild))...
%         +(a^3+a*b^2)*3/4*G*omega^alpha*sin(alpha_tild)...
%         -(omega^alpha*sin(alpha_tild));
% end