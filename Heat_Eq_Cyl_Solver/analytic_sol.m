
%-- Compute analytical solution to Heat Eq. in Cyl Coord -- 

%-------------- Heat Eq in Cylindrical Coord --------------
%----------------- Author: Owen Morehead ------------------
%------------------ Date: June 7, 2022 --------------------

function U_anal = analytic_sol(r_array,t_array,r1,r2,nu)

%-- Input:
%        r_array: column vector of spacial domain
%        t_array: array of t values at which solution will be computed
%        r1, r2: boudnary points in spacial domain
%        nu: parameter value
%
%-- Output:
%         U_anal: size(length(r_array),length(t_array))
%                 matrix of solution column vectors at each timestep
%
%--------------------------------------------------------------------------


options = optimset('Display','off');

%-- functions to be used --

Rhat = @(r,beta) ((besselj(0,r*beta).*bessely(0,r2*beta) - ...
                 besselj(0,r2*beta).*bessely(0,r*beta))./...
            ((2/(pi^2*beta^2))*((besselj(0,beta*r1).^2 - besselj(0,beta*r2).^2)./...
            besselj(0,beta*r1).^2)).^(1/2));

u0 = @(r) 10*(r-1).*(4-r).^2.*exp(-r);

int_inside = @(r,beta) r.*((besselj(0,r*beta).*bessely(0,r2*beta) - ...
                 besselj(0,r2*beta).*bessely(0,r*beta))./...
            ((2/(pi^2*beta^2))*((besselj(0,beta*r1).^2 - besselj(0,beta*r2).^2)./...
            besselj(0,beta*r1).^2)).^(1/2)).*(10*(r-r1).*(r2-r).^2.*exp(-r));
        
[x,w]=lgwt(500,r1,r2); % compute the Legendre-Gauss nodes and weights

%--------------------------------------------------------------------------
        
%-- solve eigenvalue eq for beta(k) --
for k = 1:50

%-- compute beta(k), solutions to nonlinear eq.. --

eigeq = @(beta) (besselj(0,r1*beta).*bessely(0,r2*beta) - ...
                 besselj(0,r2*beta).*bessely(0,r1*beta));
             
beta0 = 1 + 1.05*(k-1);

betas(k) = fsolve(eigeq,beta0,options);

%-----------------------------------------------------

%-- compute integrals in analytical solution --

%ints_lgwt(k) = quadl(@(r)int_inside(r, betas(k)),r1,r2);

%-- evaluating integral using lgwt.m function, Legendre-Gauss Quadrature --
int_inside_lgwt = int_inside(x,betas(k));

ints_lgwt(k) = w'*int_inside_lgwt;

end


indexer = 1;

for t = t_array

for k = 1:50
    
U_term(:,k) = exp(-betas(k)^2*nu*t).*Rhat(r_array,betas(k)).*ints_lgwt(k);

end

U_anal(:,indexer) = sum(U_term,2);

indexer = indexer + 1;

end


end

            