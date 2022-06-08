% -- 2nd ORDER FINITE DIFFERENCES METHOD TO SOLVE NONLINEAR PDE --
% -- code written by Owen Morehead --
% --     date: May 14, 2022      --

% Computing numerical soln of PDE using -> 2nd Order Finite Differences.  
% This function turns PDE into system of ODE's which is then solved using
% any desired time stepping scheme. 
% Output of this function will be the RHS of semi-discrete form

% KURAMOTO-SIVASHINSKY NONLINEAR PDE: 
% du/dt + du/dx + d''u/du'' + d''''u/du'''' = 0

% I.C: u(x,t=0) = u0 (input parameter)
% B.C: Periodic

function rhs = FDM2_KS_nonlinear_PDE(ux, t, h)

% -- Input:
%         column vector: ux
%         timestep: t (not used for this problem)
%         spacial step size: h
% -- Output:
%         rhs: RHS of semi-discrete form of PDE, to be discretized in
%              time with any time-stepping scheme.
%--------------------------------------------------------------

%-- computing first, nonlinear term in semi-discrete eq --

A1 = diag(-1*ones(length(ux(1:end-1)),1),-1) + ...
     diag(zeros(length(ux(1:end)),1),0) + ...
     diag(1*ones(length(ux(1:end-1)),1),1);
A1(end,1) = 1;  A1(1,end) = -1;

A1_mult_ux = A1*ux;

nonlinear_term = (-1/(2*h))*ux.*(A1_mult_ux);

%--------------------------------------------
%-- computing 2nd term in semi-discrete eq --

A2 = diag(1*ones(length(ux(1:end-1)),1),-1) + ...
     diag(-2*ones(length(ux(1:end)),1),0) + ...
     diag(1*ones(length(ux(1:end-1)),1),1);
A2(end,1) = 1;  A2(1,end) = 1;

term_2 = (-1/(h^2))*(A2*ux);

%--------------------------------------------
%-- computing 3rd term in semi-discrete eq --

A3 = diag(1*ones(length(ux(1:end-2)),1),-2) + ...
     diag(-4*ones(length(ux(1:end-1)),1),-1) + ...
     diag(6*ones(length(ux(1:end)),1),0) + ...  
     diag(-4*ones(length(ux(1:end-1)),1),1) + ...
     diag(1*ones(length(ux(1:end-2)),1),2);
 
 A3(end,1) = -4; A3(end-1,1) = 1; A3(end,2) = 1;
 A3(1,end) = -4; A3(1,end-1) = 1; A3(2,end) = 1;
 
 term_3 = (-1/(h^4))*(A3*ux);
 
 %--------------------------------------------
 %-- combined RHS of semi-discrete system --
 
 rhs = nonlinear_term + term_2 + term_3;
    

end
