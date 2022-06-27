function [y,t] = AB3_Method(f,y0,dt,nsteps,iosteps)

% This function computes the numerical solution to 
% the linear, first order ODE system:
% 
%          dy/dt = f(y,t) y(0)=y0
% 
% using the 3rd order explicit Adam-Bashforth (AB3) method --
%
% Input: f  -> function handle representing f(y,t)
%        y0 -> column vector of inintial conditionn
%        dt -> delta t, constant timestep size
%    nsteps -> total number of time steps ( T=dt*nsteps )
%   iosteps -> input/output steps (one time and solution pair is
%              saved into the output matrix every iosteps
%              steps.
%
% Output: y -> y -> matrix that collects the discrete solution values columnwise
%         t -> row vector of time values corresponding to solution values

t=0;
y=y0;

%-- Startup with RK3 method for first two iterations --
ys1 = RK3_Method(f,y0,dt,1,1);
y1 = ys1(:,2);

ys2 = RK3_Method(f,y1,dt,1,1);
y2 = ys2(:,2);

for ii=3:nsteps
    ts = (ii-1)*dt;
    
    y3 = y2 + dt/12*(23*f(y2,ts)-16*f(y1,(ii-2)*dt)+5*f(y0,(ii-3)*dt));
    
    if mod(ii,iosteps) == 0  %only save data when this modulus after division equals zero
        y =[y y3];           %add value y1 to array of y values
        t = [t ts+dt];       %add timestep ts+dt to array of t values
    end
    y0 = y1;  %replace previous discrete solutions at previous step with solution at next step
    y1 = y2;  
    y2 = y3;
   
end

    
end



