% -- 2nd ORDER FINITE DIFFERENCES METHOD TO SOLVE 1D HEAT EQ PDE --
% -- code written by Owen Morehead --
% --     date: May 14, 2022      --

% Computing numerical soln of 1D PDE using -> 2nd Order Finite Differences.  
% This function returns the semi-discrete form of the PDE. This IVP for a
% linear ODE can then be solved using any time-stepping scheme for IVPs.
% Output of this function will be the RHS of semi-discrete form equation.

% 1D HEAT EQUATION PDE: 
% du/dt = alpha(d''u/du'') + q(x)
% which is solved between spacial endpoints 0 and L

% I.C: u(x,t=0) = u0
% B.C: u(0,t) = g0(t) 
%      u(L,t) = gL(t)

function rhs = FDM2_1D_HEATEQ_PDE(ux, t, q, g0, gL, h, alpha)

% -- Input:
%         column vector: ux
%         boundary x values: a, b
%         # of nodes to compute discrete solution: N
% -- Output:
%         rhs: RHS of semi-discrete form of PDE, to be discretized in
%              time with any time-stepping scheme.

%h = (b-a)/(N+1); %step size, evenly spaced grid. N+2 total points
%N1 = N;
%xv = a + [1:N1]'*h;	% internal grid points. xv is a column vector.

qv = q(xv);

%---- build h(t) vector ---- in this case hv = 0
%{
hv = qv;
hv(1) = qv(1) + alpha*g0(t)/h^2;
hv(N1) = qv(N1) + alpha*gL(t)/h^2;
%}
hv = 0;  %no q(x) in PDE for this problem


%---- build differentiation matrix ----
A = diag(-2*ones(length(ux(1:end)),1),0) + ...
    diag(1*ones(length(ux(2:end)),1),-1) + ...
    diag(1*ones(length(ux(2:end)),1),1);

A = 1/(h^2)*A;


rhs = alpha*A*ux + hv;


end





