%-- Solving PDE in 2D space using Finite Differences --
%--------- and Method of Charactaristics --------------
%----------------- Driver Program ---------------------

%-------------- Author: Owen Morehead -----------------
%--------------- Date: May 25, 2022 -------------------

clear all; close all; clc;

Ns = [40,60,100,140];

indexer = 1;

for N = Ns
    

i = 1:N; j = i;

xs = 2*pi.*i/N; ys = 2*pi.*j/N;
xs = xs'; ys = ys';

[X,Y] = meshgrid(xs,ys);

h = 2*pi/N;

%----------------------------------------------------------------------
%--------------------------- FD Method --------------------------------

f1_funct = @(x,y) sin(x).*sin(y);
f2_funct = @(x,y) 1 - exp(sin(x+y));

u0_funct = @(x,y) 1/(2*pi^2)*sin(x+y).^2; %initial condition at t = 0

dt = 0.0005; 
nsteps = 1.0/dt;
iosteps = 500;

[U_fd,T_fd] = d2_fd_pde_funct(X,Y,f1_funct,f2_funct,u0_funct,h,N,dt,nsteps,iosteps);



%----------------------------------------------------------------------
%---------------- Method of Charactaristics ---------------------------


% trace charactaristic curve backwards in time to t = 0
% to give us the spacial starting point at t = 0.

f_back = @(y,t) [-sin(y(1)).*sin(y(2));
                       -(1 - exp(sin(y(2) + y(1))))];

y0_back = @(y_init) [y_init(1) ; 
                           y_init(2)];  %unknown initial condition

                       
% starting point found -> [U_back(1,end); U_back(2,end)]
% integrate u forward in time along ch curve with following ODE system

f_forw = @(y,t) [sin(y(1)).*sin(y(2));
                 (1 - exp(sin(y(2) + y(1))));
                 -(cos(y(1)).*sin(y(2)) - cos(y(1) + y(2)).*exp(sin(y(2) + y(1)))).*y(3)];
             
y0_forw = @(y_init) [y_init(1);
                           y_init(2);
                           (1/(2*pi^2))*sin(y_init(1) + y_init(2)).^2];
                       

u0_funct = @(x,y) 1/(2*pi^2)*sin(x+y).^2; %initial condition at t = 0

dt_mc = 0.01;

times = [0,0.25,0.5,0.75,1];

[U_mc,T_mc] = d2_mc_pde_funct(xs,ys,f_back,y0_back,...
    f_forw,y0_forw,u0_funct,dt_mc,times);


%----------------------------------------------------------------------

%-- calculate and plot mean squared error between numerical solutions --

%-- using quadrature rule to approx mean squared error --

mean_sq_error(indexer) = (4*pi^2/N^2)*sum(sum((U_fd(:,:,end) - U_mc(:,:,end)).^2));

indexer = indexer + 1;

end


%-- plot mean sq error vs. N  at time t = 1.0 --

figure();
plot(Ns,mean_sq_error,'linewidth',2); grid on;
xlabel('Number of Grid Points: $N$','fontsize',18,'interpreter','latex');
ylabel('Mean Squared Error','fontsize',18,'interpreter','latex');
title({sprintf('Mean Squared Error Between Numerical Solutions:'),...
       sprintf('2nd order Finite Differences and Method of Charactaristics'),...
       sprintf('At Time: t = 1.0')},...
               'interpreter','latex','fontsize',18);
 