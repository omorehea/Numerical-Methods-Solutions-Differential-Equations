% -- code written by Owen Morehead --
% --     date: May 14, 2022      --

% -- DRIVER PROGRAM --
% -- 2nd ORDER FINITE DIFFERENCES METHOD AND TIME STEPPING (AB2) TO SOLVE: 
%    initial-boundary 4th-order nonlinear PDE --

% -- PDE: Kuramoto-Sivashinsky equation --
% du/dt + du/dx + d''u/du'' + d''''u/du'''' = 0

% I.C: u(x,t=0) = u0
% B.C: Periodic

close all; clear all; clc

% -- boundaries in space domain --
a = -30; b = 30;

N = 200;   %total number of spacial points

h = 2*b/N; %dx

j = 0:N-1; 
xj = j*h + a; xj = xj'; %x array, x_0 to x_N_1
                        %does not include final boundary point
                        %since x_0 = x_N

y0 = exp(-1*(xj.^(2))); %u(x,t=0)
u0 = y0;

dt = .0005;    %timestep
T = 30;        %final time
nsteps = T/dt; %total number of timesteps

tt = 0:T/nsteps:T;

iosteps = 200; 

%{
%--------------- Euler Foward Time-stepping ---------------

t=0;
y=y0;
%y2s = y0;

for ii=1:nsteps

    ts = (ii-1)*dt;

    y1 = y0 + dt*FDM2_KS_nonlinear_PDE(y0,ts,h);
    %y2 = y0 + dt*FDM2_KS_nonlinear_PDE(y0,ts,h);
  
    if mod(ii,iosteps) == 0  %only save data when this modulus after division equals zero
        y =[y y1];       %add value y1 to array of y values
        t = [t ts+dt];   %add timestep ts+dt to array of t values
        %y2s = [y2s y2];
    end

    y0=y1; %replace discrete solution at previous step with new solution
   
end
%}

 
 
%{
%------------- AB2 Time-stepping -------------

%-- Startup with RK2 METHOD for first iteration --
t=0;
y=y0;

for ii=1:1
    
    ts = (ii-1)*dt;
    
    K1 = FDM2_KS_nonlinear_PDE(y0,ts,h);

	K2 = FDM2_KS_nonlinear_PDE(y0 + dt*K1, ts + dt,h);
    
    y1 = y0 + (dt/2)*(K1+K2);
    
  
    %if mod(ii,iosteps) == 0  %only save data when this modulus after division equals zero
        y =[y y1];      %add value y1 to array of y values
        t = [t ts+dt];  %add timestep ts+dt to array of t values
    %end
    
    y0=y1; %replace discrete solution at previous step with new solution
   
end

%-- AB2 METHOD --

y1 = y(:,2); %use 2nd solution calculated from RK2 to start AB2
y0 = y(:,1); %use 1st initial condition solution in AB2

for ii=2:nsteps
    ts = (ii-1)*dt;
    
    y2 = y1 + (dt/2)*(3*FDM2_KS_nonlinear_PDE(y1,ts,h) - ...
                      FDM2_KS_nonlinear_PDE(y0,(ii-2)*dt,h));
    
    if mod(ii,iosteps) == 0  %only save data when this modulus after division equals zero
        y =[y y2];           %add value y1 to array of y values
        t = [t ts+dt];       %add timestep ts+dt to array of t values
    end
    y0 = y1;  %replace previous discrete solutions at previous step with solution at next step
    y1 = y2;  

end

%}



%-- can acomplish AB2 timestepping using previously defined funciton --

[y,t] = AB2_Method(@(ux,t) FDM2_KS_loops(ux,t,h),y0,dt,nsteps,iosteps);




%-------------------------------------------------------------------------

y = [y;y(1,:)]; %periodic boundary.
                %append initial condition as final condition.
                
xj = [xj;(xj(end,1)+h)]'; %append periodic boudnary to x array also.

[T,X] = meshgrid(t,xj);

s = surf(X,T,y); %shading interp; view(0,90);
colorbar;
set(s,'LineStyle','none')

xlabel('x','fontsize',18,'interpreter','latex');
ylabel('t','fontsize',18,'interpreter','latex');
zlabel('u(x,t)','fontsize',18,'interpreter','latex');
title({sprintf('Numerical Solution to K-S IBVP'),sprintf('$U(x,0) = e^{-x^2}$')},...
      'fontsize',19,'interpreter','latex');


