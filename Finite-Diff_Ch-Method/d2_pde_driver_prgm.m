%-- Solving PDE in 2D space using Finite Differences --
%--------- and Method of Charactaristics --------------
%----------------- Driver Program ---------------------

%-------------- Author: Owen Morehead -----------------
%--------------- Date: May 25, 2022 -------------------

clear all; close all; clc;

%-- set parameter values --

N = 80;

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
nsteps = 2.0/dt;
iosteps = 500;

[U_fd,T_fd] = d2_fd_pde_funct(X,Y,f1_funct,f2_funct,u0_funct,h,N,dt,nsteps,iosteps);

%-- solution surface plots --
figure();

for t_i = 1:(max(size(T_fd))-3)
    %j = find(T_fd == t_i);
    subplot(3,2,t_i);
    contourf(X,Y,U_fd(:,:,t_i),20); colorbar;
    %s = surf(X,Y,U(:,:,j)); colorbar;
    %set(s,'LineStyle','none');
    
    xlabel('x','fontsize',18,'interpreter','latex');
    ylabel('y','fontsize',18,'interpreter','latex');
    title(sprintf('t = %.2f',T_fd(t_i)),'fontsize',18,'interpreter','latex');
    
    %j = j + 1;
    
end
sgtitle('Solution From Finite Differences','fontsize',20,'interpreter','latex');

%-- compute integral of F.D numerical soln vs. time --
%-- quadrature rule to represent integral --

k = 1;
for t_i = T_fd(1:2:end)
    index = find(T_fd == t_i);
    
    int_u(k) = (4*pi^2/N^2)*sum(sum(U_fd(:,:,index)));

    k = k + 1;
end

%-- plot integral results for multiple time values --
figure();
plot(T_fd(1:2:end),int_u); grid on;
xlabel('Time','fontsize',18,'interpreter','latex');
ylabel('Integral of F.D Numerical Solution','fontsize',18,'interpreter','latex');



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

times = [0,0.25,0.5,0.75,1,1.25,1.5,1.75,2.0];

[U_mc,T_mc] = d2_mc_pde_funct(xs,ys,f_back,y0_back,...
    f_forw,y0_forw,u0_funct,dt_mc,times);


 %-- plotting solution at various times --
 figure();
 for t_i = 1:(max(size(T_fd))-3)
    subplot(3,2,t_i);
    contourf(X,Y,U_mc(:,:,t_i),20); colorbar;
    xlabel('x','fontsize',18,'interpreter','latex');
    ylabel('y','fontsize',18,'interpreter','latex');
    title(sprintf('t = %.2f',T_mc(t_i)),'fontsize',18,'interpreter','latex');
 end
 sgtitle('Solution From Method of Charactaristics','fontsize',20,'interpreter','latex');


%----------------------------------------------------------------------

%-- calculate and plot max pointwise error vs. time --

for t_i = 2:max(size(times))
    
    diff_mat = abs(U_fd(:,:,t_i) - U_mc(:,:,t_i));
    error_max(t_i-1) = max(diff_mat,[],'all');
 
end

%-- plot --

figure();
semilogy(times(2:end),error_max); grid on;
xlabel('t','fontsize',18,'interpreter','latex');
ylabel('Max Pointwise Error','fontsize',18,'interpreter','latex');
title({sprintf('Max Pointwise Error Between Numerical Solutions:'),...
       sprintf('2nd Order Finite Differences and Method of Charactaristics')},...
       'interpreter','latex','fontsize',18);
 