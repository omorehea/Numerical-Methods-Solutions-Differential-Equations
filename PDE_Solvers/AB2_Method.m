function [y,t] = AB2_Method(f,y0,dt,nsteps,iosteps)

% This function computes the numerical solution to 
% the linear, first order ODE system:
% 
%          dy/dt = f(y,t) y(0)=y0
% 
% using the 2nr order explicit Adam-Bashforth (AB2) method --
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

%-- Startup with RK2 method for first iterations --
ys1 = RK2_Method(f,y0,dt,1,1);
y1 = ys1(:,2);

for ii=2:nsteps
    ts = (ii-1)*dt;
    
    y2 = y1 + (dt/2)*(3*f(y1,ts) - f(y0,(ii-2)*dt));
    
    if mod(ii,iosteps) == 0  %only save data when this modulus after division equals zero
        y =[y y2];           %add value y1 to array of y values
        t = [t ts+dt];       %add timestep ts+dt to array of t values
    end
    y0 = y1;  %replace previous discrete solutions at previous step with solution at next step
    y1 = y2;  

   
end

    
end



