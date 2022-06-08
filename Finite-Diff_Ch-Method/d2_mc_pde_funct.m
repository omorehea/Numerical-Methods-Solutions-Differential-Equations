%-- Solving PDE in 2D space using Method of Charactaristics --
%----------------- Author: Owen Morehead ---------------------
%------------------ Date: May 25, 2022 -----------------------

%clear all; close all; clc;

function [U_mc,T_mc] = d2_mc_pde_funct(xs,ys,f_back,y0_back,...
    f_forw,y0_forw,u0_funct,dt,times)

[X,Y] = meshgrid(xs,ys); %form 2d spacial grid


u0 = u0_funct(X,Y);

U_mc(:,:,1) = u0;


%-- loop over time --

for ti = 2:size(times,2)    
                     
 %   dt = 0.01;

    nsteps = times(ti)/dt;

    iosteps  = nsteps; %only need to record the final time from RK4
   
    %-- loop over all spacial domain --
    for yi = 1:size(ys)
        for xi = 1:size(xs)

            xy = ([xs(xi);ys(yi)]);

            sol = ch_method(f_back,xy,f_forw,y0_forw,dt,nsteps,iosteps);
        
            U_mc(xi,yi,ti) = sol;
        
        end
    end
    
   
    
end 

T_mc = times;

end

    
    
%--------------------------------------------------------------------------

function ch_soln = ch_method(f_back_fun, y0_back, f_forw_fun, y0_forw, dt, nsteps, iosteps)

iosteps_back  = nsteps; %only need the solution at the 'final time', t = 0

%-- backwards integration --
[U_B,T_B] = RK4_Method(f_back_fun,y0_back,dt,nsteps,iosteps_back);  
             
%-- fowards integration to get solution at time T_F --
[U_F,T_F] = RK4_Method(f_forw_fun,y0_forw([U_B(1,end),U_B(2,end)]),dt,nsteps,iosteps); 

ch_soln = U_F(end,2);

end


