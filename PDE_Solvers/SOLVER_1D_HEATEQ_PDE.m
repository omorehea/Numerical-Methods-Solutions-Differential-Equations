% -- DRIVER PROGRAM TO CALL:
% -- 2nd ORDER FINITE DIFFERENCES METHOD and
% -- GAUSS-CHEBYSHEV-LOBATTO COLLOCATION METOD
% -- TO SOLVE 1D HEAT EQUATION --


% -- code written by Owen Morehead --
% --     date: May 14, 2022      --

% IBVP: du/dt = alpha(d''u/du'') + q(x)

% which is solved between spacial endpoints a = 0 and b = L

% I.C: u(x,t=0) = u0
% B.C: u(0,t) = g0(t) 
%      u(L,t) = gL(t)

close all; clear all; clc

%define functions to be used
%syms xsym;

q =  @(x) 0*x;
g0 = @(t) 2*(t./t);
gL = @(t) 4*(t./t);

a = -1; b = 1;        %initial and final spacial boundaries
a_val = 2; b_val = 4; %constant solution value at spacial boundaries
alpha = 1;

Nvals = [25,50,100,200]; %desired number of grid points to use
%Nvals = Nvals - 2; %want N total gridpoints including ends
                   %which means subtract 2 points since algorithm uses N+2 points

I = 1; %indexing variable for appending toerror array.

for N = Nvals
%-----------------------------------------------------------------
%----- Parameters used for F-D 2 Method ------
    
h = (b-a)/(N+1); %step size, evenly spaced grid. N+2 total points
N1 = N;
xv = a + [1:N1]'*h;	% internal grid points. xv is a column vector.

%----- build differentiation matrix and RHS of 
%--    semi discrete equation : 
%--    du/dt = diff_mat*(u) + h_mat

diff_mat = diag(-2*ones(length(xv(1:end)),1),0) + ...
    diag(1*ones(length(xv(2:end)),1),-1) + ...
    diag(1*ones(length(xv(2:end)),1),1);

diff_mat = 1/(h^2)*diff_mat;

h_mat = zeros(length(xv),1);
h_mat(1) = a_val/h^2;
h_mat(end) = b_val/h^2;

y0 = (3+xv) + 5*(1-xv.^2).^2; 
u0 = [a_val;y0;b_val];

%-----------------------------------------------------------------


%----- Compute Numerical Solution using G-C-L Collocation Method ------

[x_gcl,D1,D2] = get_GCL_points_and_D_matrices(a,b,N);

y0_gcl = (3+x_gcl) + 5*(1-x_gcl.^2).^2; 

%-----------------------------------------------------------------

%-- parameters to use for time-stepping

dt = 0.0001;
T = 2;
tt = 0:dt:T;

nsteps = T/dt;
iosteps = 10;

%{
%-- Euler Foward Time-stepping --

t=0;
y=y0;

for ii=1:nsteps
    
    ts = (ii-1)*dt;
    y1 = y0 + dt*FDM2_1D_HEATEQ_PDE(y0,ts,q,g0,gL,a,b,N,alpha);
     
    %K1 = FDM2_1D_HEATEQ_PDE(y0,ts,q,g0,gL,a,b,N,alpha);
	%K2 = FDM2_1D_HEATEQ_PDE(y0 + dt*K1, ts + dt,q,g0,gL,a,b,N,alpha);
    
    %y1 = y0 + (dt/2)*(K1+K2);
    
  
    if mod(ii,iosteps) == 0  %only save data when this modulus after division equals zero
        y =[y y1];       %add value y1 to array of y values
        t = [t ts+dt];   %add timestep ts+dt to array of t values
    end
    
    y0=y1; %replace discrete solution at previous step with new solution
   
end
%}

%-----------------------------------------------------------------

%-------- Crank Nicolson Time-Stepping ---------
%-------------- FOR F.D 2 METHOD ---------------
t=0;
y=y0;

for ii = 1:nsteps
    
    ts = (ii-1)*dt;
    y1 = (eye(length(y0)) - (alpha*dt/2)*diff_mat)\...
         ((eye(length(y0)) + (alpha*dt/2)*diff_mat)*y0 + dt*h_mat);

    if mod(ii,iosteps) == 0  %only save data when this modulus after division equals zero
        y =[y y1];       %add value y1 to array of y values
        t = [t ts+dt];   %add timestep ts+dt to array of t values
    end
    
    y0 = y1;
    
end

%append initial and final condition values
y = [a_val*ones(1,length(t)); y; b_val*ones(1,length(t))]; 

%-----------------------------------------------------------------

%-------- Crank Nicolson Time-Stepping ---------
%-------------- FOR G.C.L METHOD ---------------
t_gcl=0;
y0_gcl_alt = y0_gcl(2:end-1); %only use the inner nodes for timestepping
y_gcl=y0_gcl_alt; %append initial solution vector to y_gcl

D2_alt = D2(2:end-1,2:end-1);  %inner differentiation matrix
eye_alt = eye(length(D2_alt)); 
ff_alt = 2*D2(2:end-1,1) + 4*D2(2:end-1,end); 

for ii = 1:nsteps
    
    ts = (ii-1)*dt;
    y1_gcl_alt = (eye_alt - (alpha*dt/2)*D2_alt)\...
         ((eye_alt + (alpha*dt/2)*D2_alt)*y0_gcl_alt + dt*ff_alt);

    if mod(ii,iosteps) == 0  %only save data when this modulus after division equals zero
        y_gcl =[y_gcl y1_gcl_alt];       %add value y1 to array of y values
        t_gcl = [t_gcl ts+dt];   %add timestep ts+dt to array of t values
    end
    
    y0_gcl_alt = y1_gcl_alt;
    
end

%append initial and final condition values
y_gcl = [a_val*ones(1,length(t)); y_gcl; b_val*ones(1,length(t))];

%-----------------------------------------------------------------

%----- Compute Analytical Solution -----

xv_all = [a;xv;b];
xv = [a;xv;b];

uu_fd2 = analytical_sol(xv_all,t,30);
uu_gcl = analytical_sol(x_gcl,t,30);

%------ Calculate max pointwise error ------

diff_vec_fd2 = abs(y(:,end)-uu_fd2(:,end));
error_max_fd2 = max(diff_vec_fd2,[],'all');

diff_vec_gcl = abs(y_gcl(:,end)-uu_gcl(:,end));
error_max_gcl = max(diff_vec_gcl,[],'all');

errors_fd2(I) = error_max_fd2;
errors_gcl(I) = error_max_gcl;

I = I + 1;

end

figure(1)
semilogy(Nvals,errors_fd2,'linewidth',2); grid on; hold on;
semilogy(Nvals,errors_gcl,'linewidth',2);
ylabel('Max Pointwise Error, $e(t)$','fontsize',20,'interpreter','latex');
xlabel('Number of Grid Points, $N$','fontsize',20,'interpreter','latex');
title({sprintf('Max Pointwise Error, $e(t)$, Between Analytical and Numerical Solutions'),...
       sprintf('$e(t) =$ max$_{i = 1,...,N}|U(x_i,t) - u_i(t)|$')},...
       'interpreter','latex','fontsize',18);
legend('2nd Order Finite Differences','GCL Collocation',...
       'fontsize',18,'interpreter','latex');


figure(2)
plot(xv_all,y(:,end),'r','DisplayName','2nd Order F.D'); grid on; hold on;
plot(x_gcl,y_gcl(:,end),'g','DisplayName','G.C.L Collocation');
plot(xv_all,uu_fd2(:,end),'b','DisplayName','Analytical');
%plot(xv,((3+xv) + 5*(1-xv.^2).^2),'r.');
xlabel('x','fontsize',18,'interpreter','latex');
ylabel('u(x,t=spec)','fontsize',18,'interpreter','latex');
legend


%-- surface plots for FD2, GCL, and Analytical solutions --

figure(3)
[TT,XX] = meshgrid(t,xv);
s = surf(XX,TT,y);
set(s,'LineStyle','none')
colorbar;
xlabel('x','fontsize',19,'interpreter','latex');
ylabel('t','fontsize',19,'interpreter','latex');
title({sprintf('U(x,t) : Solution to 1D Heat Equation'),...
      sprintf('$U(x,0) = (3+x) + 5(1-x^2)^2$')},...
      'fontsize',19,'interpreter','latex');


%{ 
figure(4)
s2 = surf(XX,TT,uu_fd2);
set(s2,'LineStyle','none')
xlabel('x','fontsize',18,'interpreter','latex');
ylabel('t','fontsize',18,'interpreter','latex');
title({sprintf('U(x,t) : Analytical Solution to 1D Heat Equation'),...
      sprintf('$U(x,0) = (3+x) + 5(1-x^2)^2$')},...
      'fontsize',16,'interpreter','latex');


[TT_g,XX_g] = meshgrid(t,x_gcl);
figure(5)
s2 = surf(XX_g,TT_g,y_gcl);
set(s2,'LineStyle','none')
xlabel('x','fontsize',18,'interpreter','latex');
ylabel('t','fontsize',18,'interpreter','latex');
title({sprintf('U(x,t) : Numerical Solution to 1D Heat Equation : Using GCL Collocaiton'),...
      sprintf('$U(x,0) = (3+x) + 5(1-x^2)^2$')},...
      'fontsize',16,'interpreter','latex');
%}

%------------------------------------------------------



function u_analytical = analytical_sol(xvals,t,k2)

term3 = (3 + xvals);

%-- for integration in solution --

%-- gll collocation method --
%[x_gll, w_gll] = gen_GLL_points_and_D_matrix(100,-1,1);
%u0 = (5*(1-((x_gll).^2)).^2);

%----------------------------

%-- symbolic variable method --
%u0sym = (5*(1-((xsym)^2))^2);

%----------------------------

%-- trapz() method using small x stepsize for more accuracy --
xmore = -1:.0001:1;
u0_xmore = (5*(1-((xmore).^2)).^2);

%----------------------------

for tt = t(1)
    allterm = 0;
    
    for k = 1:k2
            
        term1 =  exp(-1/4*k^2*pi^2*tt)*...
                        sin(k*pi*(xvals+1)/2);
        
        %gll_int = sum(w_gll.*(u0.*sin(k*pi*(x_gll+1)/2)));
        
        term2 = trapz(xmore,u0_xmore.*sin(k*pi*((xmore+1)/2)));
        
        %term2sym = int(((5*(1-((xsym).^2)).^2)*(sin(3*pi*((xsym+1)/2)))),-1,1)
        %term2sym = double(int((u0sym*sin(k*pi*((xsym+1)/2))),-1,1));
                         
        allterm = allterm + (term1*term2);
        
    end
    
    allterm_total = (allterm) + term3;
    
    u_analytical = allterm_total;

end


for tt = t(2:end)
    allterm = 0; allterm_total = 0;
    
    for k = 1:k2
            
        term1 =  exp(-1/4*k^2*pi^2*tt)*...
                        sin(k*pi*(xvals+1)/2);
        
        %gll_int = sum(w_gll.*(sin(k*pi*(x_gll+1)/2)));            
        
        term2 = trapz(xmore,u0_xmore.*sin(k*pi*((xmore+1)/2)));                         

        %term2sym = double(int((u0sym*sin(k*pi*((xsym+1)/2))),-1,1));
        
        allterm = allterm + (term1*term2);
        
    end
    
    allterm_total = (allterm) + term3;
    
    u_analytical = [u_analytical allterm_total];

end

end
                   