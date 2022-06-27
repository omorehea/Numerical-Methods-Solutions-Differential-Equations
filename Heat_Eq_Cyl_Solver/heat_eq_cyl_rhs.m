
%... Function to build RHS of semi-discrete form of PDE ... 

%---------- Heat Eq in Cylindrical Coord --------------
%-------------- Author: Owen Morehead -----------------
%--------------- Date: June 7, 2022 -------------------

%-- IBVP --
% dU/dt = nu(d^2U/dr^2 + 1/r(dU/dr))    nu = 1/2
% I.C -> U(r,0) = 10(r-1)(4-r)^2(e^-r)
% B.C -> U(r_1,t) = 0, U(r_2,t) = 0     r in [r_1 = 1, r_2 = 4]

function rhs = heat_eq_cyl_rhs(u,t,M_mat)

%-- Input:
%        u: column vector of solution at certain timestep
%           [u_1,...,u_N]^T. u_0 = u_N+1 = 0, not included in vector
%
%        M_mat: RHS differentiation matrix which multiplies u
%
%-- Output:
%         rhs: RHS of semi-discrete form of PDE, to be discretized
%              in time with any time-stepping scheme
%
%--------------------------------------------------------------------------

%-- build RHS --

rhs = M_mat*u;


end 

