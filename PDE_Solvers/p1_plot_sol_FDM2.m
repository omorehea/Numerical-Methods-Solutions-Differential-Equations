% -- DRIVER PROGRAM TO CALL:
% -- 2nd ORDER FINITE DIFFERENCES METHOD TO SOLVE LINEAR BVP --

% -- code written by Owen Morehead --
% --     date: May 14, 2022      --

% calls the 2-point BVP solver: "solve_FDM2.m" to solve:
% two-point BVP: -s(x)u'' - p(x)u' + q(x)u = r(x),
% B.C (1st boundary Neuman B.C): u'(a) = alf, u(b) = beta

close all; clear all; clc

%-- define functions used in BFP
s = @(x) (2 - x.^2);  %in assignment, labeled a(x)
p = @(x) (-2*x);      %in assignment, labeled a'(x)
q = @(x) (sin(pi*x)); %in assignment, labeled b(x)
r = @(x) (exp(-x));   %in assignment, labeled f(x), RHS to BVP 

a = -1; b = 1;
alf = 0; beta = 2;

Nvals = [50, 100, 500]; %desired number of grid points to use
Nvals = Nvals - 2;      %want N total gridpoints including ends
                        %which means subtract 2 points since algorithm uses N+2 points

markS = {'d','x','*'};                  
                   
figure(1);  
    I = 1;
for NN = Nvals
    [x,u] = p1_semidiscrete_FDM2(s,p,q,r,a,b,alf,beta,NN);
    
    plot(x,u,'markersize',2,'marker',markS{I}); grid on; hold on;
    
    I = I + 1;
    
end

xlabel('x','fontsize',18,'interpreter','latex');
ylabel('u(x)','fontsize',18,'interpreter','latex');
title('Numerical Solution Using 2nd Order Finite Differences','fontsize', 17,...
      'interpreter','latex');
  
for Nlab = Nvals
legendinfo{find(Nvals == Nlab)} = (sprintf('N = %.0f',Nlab + 2));
end

leg = legend(legendinfo,'fontsize',18,'interpreter','latex');
title(leg,'Total Number of Grid Points Used','interpreter','latex','fontsize',18);

