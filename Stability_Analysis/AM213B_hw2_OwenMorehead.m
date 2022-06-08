%% AM 213B: HW 2
% -- code written by Owen Morehead --
% --     date: April 25, 2022      --

%% Problem 1
%------ plot region of stability of AB3 Method ------
clear all; close all;

A = [0,10,-10;
    -100,-1,0;
     0,10,-100];
 
%A = [-1,3,-5,7;
%     0,-2,4,-6;
%     0,0,-4,6;
%     0,0,0,-16];
 
y0 = [10;10;10];

%y0 = [1;0;1;0];

theta = 0:.01:2*pi;

ewsA = eig(A);

%-- rho(exp(i*theta))/sigma(exp(i*theta)) --
z = (exp(3*1i*theta) - exp(2*1i*theta))./(23/12*exp(2*1i*theta) + 5/12 - 4/3*exp(1i*theta));

%-- Finding value of dt* that corresponds to the purely real eigenvalue ewsA(3) --
z_lambda3 = (exp(3*1i*pi) - exp(2*1i*pi))/(23/12*exp(2*1i*pi) + 5/12 - 4/3*exp(1i*pi));
dt3 = (z_lambda3 / ewsA(end));
%dt3 -- Value of dt* that corresponds to rescaling the negative real eigenvalue
%       such that it lies on the boundary of the absolute stability region.
%       This also brings the two complex eigenvalue inside the region.

%----------------------------------------------

%-- Finding value of dt* that corresponds to symmetric, complex eigenvalues ewsA(1) and ewsA(2) --
%-- Using fsolve to find where imag(rho(exp(i*theta))/(sigma(exp(i*theta))*ewsA(1))) = 0
z1_roots = fsolve(@(theta) imag((exp(3*1i*theta) - exp(2*1i*theta))./...
        (ewsA(1)*(23/12*exp(2*1i*theta) + 5/12 - 4/3*exp(1i*theta)))),[0 (pi)]);
    
z_lambda1 = (exp(3*1i*z1_roots(2)) - exp(2*1i*z1_roots(2)))/(23/12*exp(2*1i*z1_roots(2)) + 5/12 - 4/3*exp(1i*z1_roots(2)));
dt1 = z_lambda1 / ewsA(1);  
%dt1 -- Value of dt* that corresponds to rescaling the complex eigenvalues
%       such that they lie on the boundary of the absolute stability region.
%       The negative real eigenvalue is still outside the region.

rescaled_ewsA = ewsA.*(dt3);

%--------------------------------------------------------------

%-- Plot region of absolute stability highlighted with color in plot --
figure(1); 
grid on; hold on;
plot(real(z),imag(z),"HandleVisibility",'off')
fill(real(z),imag(z),'r','FaceAlpha',0.3)
%plot(real(ewsA),imag(ewsA),'r.','markersize',11)
xlabel("Re($z$)",'fontsize',16,'interpreter','latex');
ylabel("Im($z$)",'fontsize',16,'interpreter','latex');
title("Region of Absolute Stability for AB3 Method",'fontsize',18,'interpreter','latex')
%legend("Region of Absolute Stability","$\lambda_i^{*}$ -- Rescaled Eigenvalues of Matrix A",'fontsize',16,'interpreter','latex')
xlim([-1 1])
ylim([-1 1])
%axis equal

%--------------------------------------------------------------

%-- Plot region of absolute stability along with rescaled eigenvalues such
%-- that all eigenvalues are on or within boundary of abs stability --
figure(2); 
grid on; hold on;
plot(real(z),imag(z))%,"HandleVisibility",'off')
plot(real(rescaled_ewsA),imag(rescaled_ewsA),'r.','markersize',11)
xlabel("Re($z$)",'fontsize',16,'interpreter','latex');
ylabel("Im($z$)",'fontsize',16,'interpreter','latex');
title("Region of Absolute Stability for AB3 Method",'fontsize',18,'interpreter','latex')
legend("Region of Absolute Stability","$\lambda_i^{*}$ -- Rescaled Eigenvalues of Matrix A",'fontsize',16,'interpreter','latex')
xlim([-1 1])
ylim([-1 1])
%axis equal

%--------------------------------------------------------------

%--- find roots of stability polynomial with input point 
%--- inside and outside absolute stability boundary ---

%--- stability polynomial ---
%z^3 - z^2(1+z*+23/12) - z*(-4/3z + 5/12)

%--- z* = -0.2 --> inside boundary:
p_inside = [1 -74/120 -8/30 10/120];
r_inside = roots(p_inside);

%--- z* = -0.6 --> outside boundary:
p_outside = [1 18/120 -24/30 1/4];
r_outside = roots(p_outside);

x = cos(theta);
y = sin(theta);

%-- Plot of unit circle along with both sets of roots --
figure(3);
plot(x,y,"HandleVisibility",'off'); grid on; hold on
plot(real(r_inside),imag(r_inside),'r.','markersize',11)
plot(real(r_outside),imag(r_outside),'b.','markersize',11)
title({'Unit Disk and Stability Polynomial Roots','Determining Region of Absolute Stability for AB3 Method'},'interpreter','latex','fontsize',20)
xlabel('Re($z$)','interpreter','latex','fontsize',16)
ylabel('Im($z$)','interpreter','latex','fontsize',16)
xlim([-4 4]); ylim([-2 2]);

lgd = legend('$\lambda_j \Delta t = -1/5$ --- Inside Stability Boundary',...
    '$\lambda_j \Delta t = -3/5$ --- Outside Stability Boundary',...
    'interpreter','latex','fontsize',15);
title(lgd,"Roots of Stability Polynomial");

%--------------------------------------------------------------

%------ plot numerical solution to specific ODE problem ------

f = @(y,t) A*y;

T = 20;
dt = [10^(-4), real(dt3) + .0001, real(dt3) - .0001];

nsteps = T./dt;
iosteps = [2,2,2];

figure(4)
%--- run AB3 function ---

I = 1; %index for subplots
for i = 1:length(dt)
    [yk,tk] = AB3_Method(f,A,y0,dt(i),nsteps(i),iosteps(i));

    subplot(3,1,I)
    plot(tk,yk,'linewidth',2); grid on;
    
    xlabel('$t$','fontsize',16,'interpreter','latex')
    ylabel('$y(t)$','fontsize',16,'interpreter','latex')
    title(sprintf('Numerical Soluiton, AB3 Method: $\\Delta t = %2.4g$',dt(i)),'interpreter','latex','fontsize',16)
    legend('$y_1$','$y_2$','$y_3$','interpreter','latex','fontsize',15);
    I = I + 1;

end


%% Problem 2
% --- BDF3 Method ---
% --- Plot of boundary of absolute stability region ---

clear all; close all;

%-- calculate roots of first ch poly --
%-- rho(z) = z^3 - (18/11)z^2 + (9/11)z - 2/11

p1 = [1 -18/11 9/11 -2/11];
roots_p1 = roots(p1) %all roots have absolute value <= 1
                     %root at z = 1 is simple.

theta = 0:.01:2*pi;

%--- z = dt*lambda_i = rho(z)/sigma(z) ---
z = (exp(3*1i*theta) - 18/11*exp(2*1i*theta) + 9/11*exp(1i*theta) - 2/11)./(6/11*exp(3*1i*theta));


%-- Plot region of absolute stability highlighted with color in plot --

figure(1); 
plot(real(z),imag(z)); grid on
xlabel("Re($z$)",'fontsize',16,'interpreter','latex');
ylabel("Im($z$)",'fontsize',16,'interpreter','latex');
title("Region of Absolute Stability for BDF3 Method",'fontsize',18,'interpreter','latex')
ylim([-5 5])
xlim([-2 7])
text(-1.5,4.2,'Region of Stability, $R$','interpreter','latex','fontsize',16)

%--------------------------------------------------------------

%--- find roots of stability polynomial with input point 
%--- inside and outside absolute stability boundary ---

%--- z* = 3.0 --> inside boundary:
p_inside = [-7/11 -18/11 9/11 -2/11];
r_inside = roots(p_inside);

%--- z* = -1.0 --> inside boundary:
p_outside = [17/11 -18/11 9/11 -2/11];
r_outside = roots(p_outside);

%--------------------------------------------------------------

x = cos(theta);
y = sin(theta);

%-- Plot of unit circle along with both sets of roots --

figure(2);
plot(x,y,"HandleVisibility",'off'); grid on; hold on
plot(real(r_inside),imag(r_inside),'r.','markersize',11)
plot(real(r_outside),imag(r_outside),'b.','markersize',11)
title({'Unit Disk and Stability Polynomial Roots','Determining Region of Absolute Stability for BDF3 Method'},'interpreter','latex','fontsize',20)
xlabel('Re($z$)','interpreter','latex','fontsize',16)
ylabel('Im($z$)','interpreter','latex','fontsize',16)
xlim([-4 4]); ylim([-2 2]);

lgd = legend('$\lambda_j \Delta t = 3.0$ --- Inside Stability Boundary',...
    '$\lambda_j \Delta t = -1.0$ --- Outside Stability Boundary',...
    'interpreter','latex','fontsize',15);
title(lgd,"Roots of Stability Polynomial");

%we see the region of abs stability is everything outside the booundary

%{
figure(4); 
grid on; hold on;
plot(real(z),imag(z),"HandleVisibility",'off')
fill(real(z),imag(z),'w')
%plot(real(rescaled_ewsA),imag(rescaled_ewsA),'r.','markersize',11)
xlabel("Re($z$)",'fontsize',16,'interpreter','latex');
ylabel("Im($z$)",'fontsize',16,'interpreter','latex');
title("Region of Absolute Stability for AB3 Method",'fontsize',18,'interpreter','latex')
ylim([-4.5 4.5])

ax = gca;
ax.Color = [.8 .2 0];
%}

%% Problem 3
% Numerical Scheme:
% --- u_k+2 - 4u_k+1 + 3u_k = -2(dt)f(u_k,t_k) ---

% --- Plot of boundary of absolute stability region ---

clear all; close all;

theta = 0:.01:2*pi;

z = (4*exp(1i*theta) - exp(2*1i*theta) - 3)/2;

%--- Plot region of absolute stability boundary ---
figure(1);
plot(real(z),imag(z)); grid on;
xlabel("Re($z$)",'fontsize',16,'interpreter','latex');
ylabel("Im($z$)",'fontsize',16,'interpreter','latex');
title({sprintf("Region of Absolute Stability"),sprintf("Explicit Linear Multistep Scheme -- HW2 P3")},'fontsize',18,'interpreter','latex')
xlim([-5 1]); ylim([-4 4]);

%--------------------------------------------------------------

%--- find roots of stability polynomial with input point 
%--- inside and outside absolute stability boundary ---

%--- z* = -2.0 --> inside boundary:
p_inside = [1 -4 -1];
r_inside = roots(p_inside);

%--- z* = -8.0 --> inside boundary:
p_outside = [1 -4 -13];
r_outside = roots(p_outside);

%--------------------------------------------------------------

x = cos(theta);
y = sin(theta);

%-- Plot of unit circle along with both sets of roots --

figure(2);
plot(x,y,"HandleVisibility",'off'); grid on; hold on
plot(real(r_inside),imag(r_inside),'r.','markersize',11)
plot(real(r_outside),imag(r_outside),'b.','markersize',11)
title({'Unit Disk and Stability Polynomial Roots','Determining Region of Absolute Stability'},'fontsize',20,'interpreter','latex')
xlabel('Re($z$)','interpreter','latex','fontsize',16)
ylabel('Im($z$)','interpreter','latex','fontsize',16)
%xlim([-4 4]); ylim([-2 2]);
axis equal

lgd = legend('$\lambda_j \Delta t = -2.0$ --- Inside Stability Boundary',...
    '$\lambda_j \Delta t = -8.0$ --- Outside Stability Boundary',...
    'interpreter','latex','fontsize',15);
title(lgd,"Roots of Stability Polynomial");

%we see this method is not absolutely stable.


%% Problem 4
% Numerical Scheme: Implicit Midpoint Method
% --- u_k+1 - u_k = (dt)f((u_k+1 + u_k)/2,t_k + dt) ---

% --- Plot of boundary of absolute stability region ---
% --- Abs stability region is everywhere where Re(z) < 0 ---

clear all; close all;

xs = [-10,0,0,-10]; ys = [-10,-10,10,10];

% --plot region of absolute stability highlighted in red --

figure(1); 
fill(xs,ys,'r','FaceAlpha',0.3); grid on
xlabel("Re($z$)",'fontsize',16,'interpreter','latex');
ylabel("Im($z$)",'fontsize',16,'interpreter','latex');
title({sprintf("Region of Absolute Stability"),sprintf("Implicit Midpoint Method")},'fontsize',18,'interpreter','latex')
xlim([-2 2]); ylim([-2 2]);
%text(-1.5,4.2,'Region of Stability, $R$','interpreter','latex','fontsize',16)


%{
%--- z* = -2.0 --> inside boundary:
p_inside = [1 -4 -1];
r_inside = roots(p_inside)

%--- z* = -8.0 --> inside boundary:
p_outside = [1 -4 -13];
r_outside = roots(p_outside)

figure(2);

x = cos(theta);
y = sin(theta);

plot(x,y,"HandleVisibility",'off'); grid on; hold on
plot(real(r_inside),imag(r_inside),'r.','markersize',11)
plot(real(r_outside),imag(r_outside),'b.','markersize',11)
title({'Unit Disk and Stability Polynomial Roots','Determining Region of Absolute Stability'},'fontsize',20,'interpreter','latex')
xlabel('Re($z$)','interpreter','latex','fontsize',16)
ylabel('Im($z$)','interpreter','latex','fontsize',16)
%xlim([-4 4]); ylim([-2 2]);
axis equal

lgd = legend('$\lambda_j \Delta t = -2.0$ --- Inside Stability Boundary',...
    '$\lambda_j \Delta t = -8.0$ --- Outside Stability Boundary',...
    'interpreter','latex','fontsize',15);
title(lgd,"Roots of Stability Polynomial");


%}









 