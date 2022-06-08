function [y,t] = RK4_Method(f,y0,dt,nsteps,iosteps)

% This function computes the numerical solution to 
% the linear, first order ODE system:
% 
%          dy/dt = f(y,t) , y(0)=y0
%
% using the explicit 4th stage Runge-Kutta (RK4) method --
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

	K2 = f(y0+ dt * 0.5 *K1, ts + 0.5* dt);

	K3 = f(y0+ dt*0.5*K2, ts + 0.5* dt);
    
    K4 = f(y0+ dt*K3, ts + dt);
    
    y1 = y0 + (dt/6)*(K1 + 2*K2 + 2*K3 + K4);
    
  
    if mod(ii,iosteps) == 0  %only save data when this modulus after division equals zero
        y =[y y1];      %add value y1 to array of y values
        t = [t ts+dt];  %add timestep ts+dt to array of t values
    end
   y0=y1; %replace discrete solution at previous step with new solution
   
end

    
end



