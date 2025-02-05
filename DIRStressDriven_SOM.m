function [t_vec, strain_vec] = DIRStressDriven_SOM(par, stress, tspan, strain_0)
% DIRStressDriven_SOM
% Single Order Material Model of a Viscoelastic Material, time-domain response
% Based on the Diffusive Representations of Fractional Integrals method.
% 
% Author: Dr. Afshin Farhadi, THWS, Germany
% Email: afshin.farhadi19@gmail.com
% Created: 1 November 2024
% 
% Inputs:
% par        ... (4-by-1)-array of parameters:
%                 alpha = par(1) \in (0,1)
%                 E0 = par(2)
%                 E1 = par(3)
%                 p1 = par(4)
% stress     ... function handle for stress(t)
% tspan      ... time span [t_0, t_1]
% strain_0   ... initial condition of strain (scalar)
% 
% Outputs:
% t_vec      ... time vector
% strain_vec ... strain solution vector corresponding to t_vec

    % Unpack parameters
    alpha = par(1);
    E0 = par(2);
    E1 = par(3);
    p1 = par(4);

    % Derived parameters
    b = E1 / p1;
    c = E0 + E1;
    d = E0 * E1 / p1;
    calpha = sin(pi * alpha) / pi;

    % Initialize time grid
    tmin = tspan(1);
    tmax = tspan(2);
    

close all; clc

NewtonFlag = 2;          % 1: Newton Raphson, 2: Modified Newton Raphson

Ndim =1;                 % dimension of the model

maxit=100;

initial_sigma = stress(tmin);           % initial condition of stress

initial_epsilon = strain_0;         % initial condition of strain

calpha=sin(pi*alpha)/pi;

Nt = 10^5;                   % Number of time grid points

%h=10^-5;

N  = 150;                               % Number of Gauss-Laguerre quadrature points

tolerance=1e-10;                        % Tolerance for the Newton's method

tmin        = 0;                        % Lower bound of the domain

tmax        = 1;                        % Upper bound of the domain

h           = (tmax-tmin)/(Nt-1);       % Average grid size

t_vec   = (tmin:h:tmax);                    % Coordinate of time grid nodes


% Nt = ceil((tmax-tmin)/h);             % Number of time grid points
% % 
% t(1) = 0;
% % 
% for i=2:Nt
% 
%     t(i) = t(i-1)+h;
% 
% end 

[Xnode,W]=Gaulagwt(N);                 % function for computing nodes and weights

%--------------------------------------------------------------------------
% We suppose that stress sigma, as a function of time, is given as follows:

sigma= stress;

%--------------------------------------------------------------------------
% Defining auxiliary terms
%-----------------------------------------
r1=zeros(Ndim,N);

r2=zeros(Ndim,N);

for j=1:Ndim
    
    for i=1:N

        r1(j,i)=(-1/(1-alpha(j))).*Xnode(i);

        r2(j,i)=(1/alpha(j)).*Xnode(i);
    end
end

for j=1:Ndim

     for l=1:N

          V_1(j,l) =  1./(1+(h*exp(r1(j,l))));

          V_2(j,l) = V_1(j,l)*h*calpha(j)*exp((1-alpha(j))*r1(j,l));

          V_3(j,l) = exp(-r2(j,l))/(h+exp(-r2(j,l)));

          V_4(j,l) =h*calpha(j)*exp(-alpha(j)*r2(j,l))/(h+exp(-r2(j,l)));

          V_5(j,l) = ((V_2(j,l)./(1-alpha(j)))+(V_4(j,l)./(alpha(j))));

          Z10(j,l) = 0;
    
          Z20(j,l) = 0;

          E10(j,l) = 0;
    
          E20(j,l) = 0;
     end
end


for l=1:N

    Xnode(l)=W(l)*exp(Xnode(l));

end 

strain_vec=zeros(Ndim,Nt);

for j=1:Ndim

    strain_vec(j,1)=initial_epsilon(j);

end
 
% 
J =zeros(Ndim,Ndim);
ynew = zeros(1,Ndim);
g1 = zeros(1,Ndim);
for k = 2:Nt                     %numerical integration as time is increased 
                                 %t vector above hold incremental time values.
    for j=1:Ndim

        strain_vec(j,k) = strain_vec(j,k-1);

        ytemp      = strain_vec(j,k); 

    end

     eps=1;

     i=1;  
           
     while eps>=tolerance 
                 
          stress   = sigma(t_vec(k));

          df1= 1;
          
          
   
         for j=1:Ndim

             int_stress = 0;

             int_strain = 0;

             for l=1:N

                  Z11(j,l) = V_1(j,l)*Z10(j,l)+V_2(j,l)*stress;

                  Z21(j,l) = V_3(j,l)*Z20(j,l)+V_4(j,l)*stress;

                  E11(j,l) = V_1(j,l)*E10(j,l)+V_2(j,l)*ytemp(j);

                  E21(j,l) = V_3(j,l)*E20(j,l)+V_4(j,l)*ytemp(j);

                  int_stress = int_stress+Xnode(l)*((1/(1-alpha(j)))*Z11(j,l)+(1/(alpha(j)))*Z21(j,l));

                  int_strain = int_strain+Xnode(l)*((1/(1-alpha(j)))*E11(j,l)+(1/(alpha(j)))*E21(j,l));

              end 
                  
              g1(j)=ytemp(j)+(d/c)*int_strain-(b/c)*int_stress-(1/c)*stress+(1/c)*initial_sigma-initial_epsilon(j);
                            
              if  i==1   || NewtonFlag  == 1      

                  for jj = 1:Ndim
                 
                      DIFR_temp = 0 ;

                      for l=1:N

                           DIFR_temp= DIFR_temp+Xnode(l)*V_5(j,l)*df1(j,jj);
                      end 

                      if j == jj

                         J(j,jj) =  1 - (d/c)*DIFR_temp;

                      else

                      J(j,jj) =   - (d/c)*DIFR_temp;

                      end
                  end
              else
                  continue
              end
          end
                                 
          if NewtonFlag == 1 

             ynew=ytemp-(g1/J);

          elseif i == 1 

                 Jinv1 = inv(J);

                 ynew = ytemp - Jinv1 * transpose(g1);
       
          else 

                 ynew = ytemp - Jinv1 * transpose(g1);
          end 

          %eps=norm(ynew-ytemp);

          eps=norm(ynew-ytemp)/norm(ytemp);      % rel Tolerance

          ytemp=ynew;

          i=i+1; 

          if  i > maxit

          disp('stop');

          error('Newton''s Method did not converge within Maxit iterations');
          end
    end
      % fprintf('Iteration Number for k= %d',k);
      % 
      % fprintf(' is equal to %d\n',i);
      
    for j=1:Ndim

        strain_vec(j,k)=ytemp(j);

        for l=1:N

            Z10(j,l) = Z11(j,l);

            Z20(j,l) = Z21(j,l);

            E10(j,l) = E11(j,l);

            E20(j,l) = E21(j,l);

         end  
     end
end
 
function [x,w] = Gaulagwt(n, alf)
% Gaulagwt.m
% This script is for computing definite integrals using Laguerre-Gauss 
% Quadrature. Computes the Laguerre-Gauss nodes and weights  on an interval
% [0, inf] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [0, inf]
% which you can evaluate at any x in [0,inf]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral using w*f;
% Input
% n   =  number of nodes for the integration
% alf = the power of x in the generalized form of gauss laguerre integral
% Output
% x   = absicsa at which the function is to be evaluated
% w   = weights corresponding to each abscisca
%
%
% Note: w is a row vector while x is a column vector. This choice is made
% so that if f is vectorized, we can compute the integration by w*f(x)
% Written by Lateef Kareem - 17/10/2018
if nargin == 1
    alf = 0;
end
MAXIT = 20;
p2 = 0;
pp = 0;
z = 0; 
eps = 3e-14;
x = zeros(n, 1); 
w = zeros(1, n);
for i = 1:n
    % Loop over the desired roots.
    if (i == 1)
        % Initial guess for the smallest root.
        z = (1.0 + alf) * (3.0 + 0.92 * alf) / (1.0 + 2.4 * n + 1.8 * alf);
    elseif (i == 2)
        %Initial guess for the second root.
        z = z + (15.0 + 6.25 * alf) / (1.0 + 0.9 * alf + 2.5 * n);
    else
        % Initial guess for the other roots.
        ai = i - 2;
        z = z + ((1.0 + 2.55 * ai) / (1.9 * ai) + 1.26 * ai * alf / ....
                (1.0 + 3.5 * ai)) * (z - x(i - 2)) / (1.0 + 0.3 * alf);
    end
    for its = 1: MAXIT
        % Refinement by Newton?s method.
        p1 = 1.0;
        p2 = 0.0;
        for j = 1:n
            % Loop up the recurrence relation to get the
            p3 = p2; % Laguerre polynomial evaluated at z.
            p2 = p1;
            p1 = ((2 * j - 1 + alf - z) * p2 - (j - 1 + alf) * p3) / j;
        end
        % p1 is now the desired Laguerre polynomial. We next compute pp, its derivative,
        % by a standard relation involving also p2, the polynomial of one lower order.
        pp = (n * p1 - (n + alf) * p2) / z;
        z1 = z;
        z = z1 - p1 / pp; % Newton?s formula.
        if (abs(z - z1) <= eps); break; end
    end
    x(i) = z; % Store the root and the weight.
    w(i) = -exp(gammaln(alf + n) - gammaln(n)) / (pp * n * p2);
end
