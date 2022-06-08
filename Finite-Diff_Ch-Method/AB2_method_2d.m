function [y,t] = AB2_method_2d(f,y0,DT,NSTEPS,IOSTEPS)

% This function computes the numerical solution to 
% the system of ordinary differential equations
% 
%          dy/dt = f(y,t) y(0)=y0
%
% Input f  -> function handle representing f(y,t)
%       y0 -> column vector of inintial conditionn
%       DT -> delta t
%   NSTEPS -> total number of time steps ( T=DT*NSTEPS)
%  IOSTEPS -> input/output steps (one time snapshots is
%             saved into the output matrix every IOSTEPS
%             steps.
%
% Output y -> matrix that collects the time snapshots of 
%             the solution columnwise
%        t -> row vector collecting the times at which the 
%             solution is saved

%-- only difference for 2D code is how the solution y(:,:,:) is built up
%-- solutions stored in a 3D matrix instead of 2D matrix

t=0;
y=y0;

% Startup with Heun method
ts=0;
y1 = y0+0.5*DT*(f(y0+DT*f(y0,ts),ts+DT)+f(y0,ts));

indexer = 2;

for ii=2:NSTEPS
    ts = (ii-1)*DT;
    
    y2 = y1 + 0.5*DT*(3*f(y1,ts)- f(y0,ts-DT));
    
    if mod(ii,IOSTEPS)==0
        y(:,:,indexer) = y2;
        t = [t ts+DT];
        
        indexer = indexer + 1;
    end
    y0 = y1; 
    y1 = y2; 
end

    
end



