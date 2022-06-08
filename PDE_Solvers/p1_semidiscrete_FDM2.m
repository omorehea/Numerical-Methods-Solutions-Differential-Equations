% -- 2nd ORDER FINITE DIFFERENCES METHOD TO SOLVE LINEAR BVP --
% -- code written by Owen Morehead --
% --     date: May 14, 2022      --

% Computing numerical soln of BVP using -> 2nd Order Finite Differences 

% two-point BVP: -s(x)u'' - p(x)u' + q(x)u = r(x),
% B.C (1st boundary Neuman B.C): u'(a) = alf, u(b) = beta

function [xb,ub] = solve_FDM2(s, p, q, r, a, b, alf, beta, N)

% -- Input:
%         functions of x: s, p, q, r
%         boundary x values: a, b
%         solution at boundarys: alf, beta
%         # of sub-intervals: N
% -- Output:
%         xb: x array with bc
%         ub: u array with bc
%

h = (b-a)/(N+1); %step size, evenly spaced grid. N+2 total points
N1 = N;
xv = a + [0:N1]'*h;	% internal grid points, including x0 = a. xv is a column vector.
sv = s(xv); % function s at internal grid points, column vector
pv = p(xv);	% function p at internal grid points, column vector
qv = q(xv);	% function q at internal grid points, column vector
rv = r(xv);	% function r at internal grid points, column vector

%-- Build vector bv
gv=rv;      % RHS of linear system, column vector
gv(1)=alf;
gv(end)=gv(end)+(sv(end)/(h^2) + pv(end)/(2*h))*beta;  

%-- Build matrix A (Derivative matrix) using diag()

A=diag((-sv(2:end)/h^2)+pv(2:end)/(2*h),-1)+diag((2*sv(1:end))/h^2 + qv(1:end),0)+...
  diag(-((sv(1:end-1)/h^2) + (pv(1:end-1)/(2*h))),1);

%-- Forward F.D used to encorporate Neuman first B.C
%-- we change the first three matrix entries in row 1 to account for this
A(1,1) = 3/(2*h); A(1,2) = -2/h; A(1,3) = 1/(2*h);

%-- Solve linear system for u
uv=A\gv;

xb=[xv; b];	% x with bc
ub=[uv; beta];	% u with bc 

end

