
% Function to build RHS differentiation matrix  for... 

%---------- Heat Eq in Cylindrical Coord --------------
%-------------- Author: Owen Morehead -----------------
%--------------- Date: June 7, 2022 -------------------

%-- IBVP --
% dU/dt = nu(d^2U/dr^2 + 1/r(dU/dr))    nu = 1/2
% I.C -> U(r,0) = 10(r-1)(4-r)^2(e^-r)
% B.C -> U(r_1,t) = 0, U(r_2,t) = 0     r in [r_1 = 1, r_2 = 4]

function M_mat = heat_eq_cyl_M_mat(u,j,dr,nu)

%-- Input:
%        u: column vector of solution at certain timestep
%        j: array of index values from 1 to n (same size as u)
%        t: array of time values (not used in this implementation)
%        dr: spacial step size in r
%        nu: parameter value
%
%-- Output:
%         M_mat: differentiation matrix in RHS of semi-discrete form of PDE
%
%--------------------------------------------------------------------------

%-- build differentiation matrix --

M2 = diag(-2*ones(length(u(1:end)),1),0) + ...
     diag(1*ones(length(u(2:end)),1),-1) + ...
     diag(1*ones(length(u(2:end)),1),1);
 
M2 = M2/dr^2;

 
M1 = diag(1./(1+(j(1:end-1))*dr),1) + ...
     diag(-1./(1+ (j(2:end))*dr),-1);

M1 = M1/(2*dr);
 
M_mat = nu*(M1 + M2);


end 

