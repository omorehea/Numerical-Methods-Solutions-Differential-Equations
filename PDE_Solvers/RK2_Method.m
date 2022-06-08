function [y,t] = RK2_Method(f,y0,dt,nsteps,iosteps)

% This function computes the numerical solution to 
% the linear, first order ODE system:
% 
%          dy/dt = f(y,t) , y(0)=y0
%
% using the explicit 2nd stage Runge-Kutta (RK2) method (Heun Method) --
%
% Input: f  -> function handle representing f(y,t)
%        y0 -> column vector of inintial condition
%        dt -> delta t, constant timestep size
%    nsteps -> total number of time steps ( T=dt*nsteps )
%   iosteps -> input/output steps (one time and solution pair is
%              saved into the output matrix every iosteps
%              steps.
%
% Output: y -> matrix that collects the discrete solution values columnwise
%         t -> row vector of time values corresponding to solution values

t=0;
y=y0;

for ii=1:nsteps
    
    ts = (ii-1)*dt;
    
    K1 = f(y0,ts);

	K2 = f(y0 + dt*K1, ts + dt);
    
    y1 = y0 + (dt/2)*(K1+K2);
    
  
    if mod(ii,iosteps) == 0  %only save data when this modulus after division equals zero
        y =[y y1];      %add value y1 to array of y values
        t = [t ts+dt];  %add timestep ts+dt to array of t values
    end
    
    y0=y1; %replace discrete solution at previous step with new solution
   
end

    
end



