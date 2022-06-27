%%
%---------- Heat Eq in Cylindrical Coord --------------
%-------------- Author: Owen Morehead -----------------
%--------------- Date: June 7, 2022 -------------------

%-- IBVP --
% dU/dt = nu(d^2U/dr^2 + 1/r(dU/dr))    nu = 1/2
% I.C -> U(r,0) = 10(r-1)(4-r)^2(e^-r)
% B.C -> U(r_1,t) = 0, U(r_2,t) = 0     r in [r_1 = 1, r_2 = 4]


clear all; close all; clc;

%-- initial parameters --

%-- these parameter values are also defined in each code section --
%-- so they can be run individually. Can also delete those and make sure --
%-- to run this first section so the parameters initialize --

r1 = 1; r2 = 4; %initial and final spacial boundaries
nu = 1/2;

u0_funct = @(r) 10*(r-1).*(4-r).^2.*exp(-r);

%--------------------------------------------------------------------------

%%
%---------------- Compute and plot analytical solution --------------------

r1 = 1; r2 = 4; %initial and final spacial boundaries
nu = 1/2;


n = 200;
j = [0:n+1];           
j_noends = j(2:end-1);
dr = (r2 - r1)/(n+1);   %spacial step size

rs_all = r1 + j*dr; rs_all = rs_all'; %array of spacial values

t_array = [0,0.5,1,2];  %array of times to evaluate at


U_anal = analytic_sol(rs_all,t_array,r1,r2,nu); %compute analytical soluiton


%figure()
%s = surf(t_array,rs_all,U_anal);
%set(s,'LineStyle','none');
%colorbar;
%xlabel('t','fontsize',18,'interpreter','latex');
%ylabel('r','fontsize',18,'interpreter','latex');


figure()

plot(rs_all,U_anal,'linewidth',1.5); grid on;
xlabel('$r$','fontsize',20,'interpreter','latex');
ylabel('$U(r,t)$','fontsize',20,'interpreter','latex');
title({('Analytical Solution to Heat Eq. in Cyl Coord.')},...
        'fontsize',20,'interpreter','latex');

for ts = t_array
legendinfo{find(ts == t_array)} = (sprintf('t = %.1f',ts));
end

legend(legendinfo,'fontsize',18,'interpreter','latex')


%--------------------------------------------------------------------------

%%
%------------------- compute spectral radius of M_mat ---------------------
%---------------------- and plot: rho(M_mat) vs. n ------------------------
%-------------------- and plot: log(dt*) vs. log(n) -----------------------

r1 = 1; r2 = 4; %initial and final spacial boundaries
nu = 1/2;

u0_funct = @(r) 10*(r-1).*(4-r).^2.*exp(-r);

k = 1:20;
count = 1;

for n = (10 + 10*(k-1)) %values of n to loop through

j = [0:n+1];

j_noends = j(2:end-1);

dr = (r2 - r1)/(n+1);

rs_all = r1 + j*dr; rs_all = rs_all';

rs_inner = rs_all(2:end-1);

u0 = u0_funct(rs_inner);  %initial condition evaluated at inner spacial nodes


M = heat_eq_cyl_M_mat(u0,j_noends,dr,nu);  %compute differentiation matrix

spec_rad = max(abs(eig(M))); %max eigenvalue of M

n_nodes(count) = n;  %record n value in array (for plotting)

spec_rad_nodes(count) = spec_rad; %record spectral radius value in array

dt_star(count) = 6/(11*spec_rad); %record dt^* value in array

count = count + 1;

end

figure()
plot(n_nodes,spec_rad_nodes,'linewidth',2); grid on;
xlabel('number of inner nodes: $n$','fontsize',18,'interpreter','latex');
ylabel('$\rho(M)$','fontsize',18,'interpreter','latex');
title('Spectral Radius of $M$ vs. Number of Inner Nodes','fontsize',18,'interpreter','latex');

figure()
loglog(n_nodes,dt_star,'linewidth',2); grid on;
xlabel('number of inner nodes: $n$','fontsize',18,'interpreter','latex');
ylabel('$\Delta t^{*}$','fontsize',18,'interpreter','latex');
title('$\Delta t^{*}$ vs. Number of Inner Nodes','fontsize',18,'interpreter','latex');


%-------------------------------------------------------------------------

%%
%--------------------- Compute numerical solution ------------------------- 
%------------- using 2nd Order Centered Finite-Differences ----------------

r1 = 1; r2 = 4; %initial and final spacial boundaries
nu = 1/2;

u0_funct = @(r) 10*(r-1).*(4-r).^2.*exp(-r);

n = 200;

j = [0:n+1];

j_noends = j(2:end-1);

dr = (r2 - r1)/(n+1);

rs_all = r1 + j*dr; rs_all = rs_all';

rs_inner = rs_all(2:end-1);

u0 = u0_funct(rs_inner);



dt = 5*10^(-5);       %timestep
t_final = 2;          %final time
nsteps = t_final/dt;  %number of time steps
iosteps = 500;


M = heat_eq_cyl_M_mat(u0,j_noends,dr,nu);  %get differentiation matrix 

[U,T] = AB3_Method(@(u,t) heat_eq_cyl_rhs(u,t,M),...  %timestep using AB3
                          u0,dt,nsteps,iosteps);

U_all = [0*ones(1,size(U,2));U;0*ones(1,size(U,2))];  %append back boundary values

[TT,RR] = meshgrid(T,rs_all);

figure()

s = surf(TT,RR,U_all);     %plot of numerical solution in space and time
set(s,'LineStyle','none');
colorbar;
xlabel('$t$','fontsize',18,'interpreter','latex');
ylabel('$r$','fontsize',18,'interpreter','latex');
zlabel('$U(r,t)$','fontsize',18,'interpreter','latex');
title({('Numerical Solution to Heat Eq. in Cyl Coord.'),...
       ('Using 2nd Order Centered Finite Differences')},...
        'fontsize',18,'interpreter','latex');
  
%--------------------------------------------------------------------------

%%
%------------ Compute and plot max pointwise error vs. time --------------- 
%--------------------- for n = [25,50,100,150] ----------------------------
    
r1 = 1; r2 = 4; %initial and final spacial boundaries
nu = 1/2;

u0_funct = @(r) 10*(r-1).*(4-r).^2.*exp(-r);

dt = 5*10^(-5);       %timestep
t_final = 2;          %final time
nsteps = t_final/dt;  %number of time steps
iosteps = 80;


%-- time dependent max point-wise error --

Ns = [25,50,100,150];

index = 1;

for n = Ns
    
j = [0:n+1];

j_noends = j(2:end-1);

dr = (r2 - r1)/(n+1);

rs_all = r1 + j*dr; rs_all = rs_all';

rs_inner = rs_all(2:end-1);

u0 = u0_funct(rs_inner);


M = heat_eq_cyl_M_mat(u0,j_noends,dr,nu);

[U,T] = AB3_Method(@(u,t) heat_eq_cyl_rhs(u,t,M),...
                          u0,dt,nsteps,iosteps);

%both numerical and analytical solutions evaluated at same spacial and time points                     
                      
U_all = [zeros(1,size(U,2));U;zeros(1,size(U,2))]; %numerical solution

U_anal = analytic_sol(rs_all,T,r1,r2,nu);          %analytical solution

U_diff = abs(U_all(2:end-1,:)-U_anal(2:end-1,:));  %difference between soluitons

pointwise_error(index,:) = max(U_diff);            %calculate error at each time value
%pointwise_error(index) = max(U_diff,[],'all');

%[j1,j2] = find(pointwise_error(index) == U_diff);

index = index + 1;

end

figure()

semilogy(T,pointwise_error,'linewidth',1.5); grid on;
xlabel('$t$','fontsize',19,'interpreter','latex');
ylabel('$e_n(t)$','fontsize',19,'interpreter','latex');
title({('Max Point-Wise Error'),('Between Numerical and Analytical Solutions')},...
        'fontsize',18,'interpreter','latex');
    

for n_val = Ns
legendinfo{find(n_val == Ns)} = (sprintf('N = %.1f',n_val));
end

legend(legendinfo,'fontsize',18,'interpreter','latex')

%--------------------------------------------------------------------------    
    
%%
%------------ Compute and plot max pointwise error at t = 2 --------------- 
%------------------------ as a function of n ------------------------------


r1 = 1; r2 = 4; %initial and final spacial boundaries
nu = 1/2;

u0_funct = @(r) 10*(r-1).*(4-r).^2.*exp(-r);

dt = 5*10^(-5);       %timestep
t_final = 2;          %final time
nsteps = t_final/dt;  %number of time steps
iosteps = 80;

%-- plot error at final time vs. n

k = 1:20;
index = 1;

Ns_ks = (10 + 10*(k-1));

for n = Ns_ks
    
    
j = [0:n+1];

j_noends = j(2:end-1);

dr = (r2 - r1)/(n+1);

rs_all = r1 + j*dr; rs_all = rs_all';

rs_inner = rs_all(2:end-1);

u0 = u0_funct(rs_inner);


M = heat_eq_cyl_M_mat(u0,j_noends,dr,nu);

[U,T] = AB3_Method(@(u,t) heat_eq_cyl_rhs(u,t,M),...
                          u0,dt,nsteps,iosteps);

%both numerical and analytical solutions evaluated at same spacial and time points                     
                      
U_all = [0*ones(1,size(U,2));U;0*ones(1,size(U,2))]; %numerical solution

U_anal = analytic_sol(rs_all,T,r1,r2,nu);            %analytical solution


U_diff = abs(U_all(2:end-1,end)-U_anal(2:end-1,end));    
    
pointwise_error_end(index) = max(U_diff,[],'all');

%U_diff(U_diff == pointwise_error(index))

index = index + 1;

end

figure()

loglog(Ns_ks,pointwise_error_end,'linewidth',1.5); grid on;
xlabel('number of inner spacial nodes: $n$','fontsize',18,'interpreter','latex');
ylabel('$e_n(2)$','fontsize',18,'interpreter','latex');
title({('Max Point-Wise Error at time $t = 2.0$'),('Between Numerical and Analytical Solutions')},...
        'fontsize',18,'interpreter','latex');
   

