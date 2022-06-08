%-- Solving PDE in 2D space using Finite Differences --
%-------------- Author: Owen Morehead -----------------
%--------------- Date: May 25, 2022 -------------------


function [U_fd,T_fd] = d2_fd_pde_funct(X,Y,f1_funct,f2_funct,u0_funct,h,N,dt,nsteps,iosteps)


f1 = f1_funct(X,Y);  %sin(X).*sin(Y);
f2 = f2_funct(X,Y);  %1 - exp(sin(X+Y));

u0 = u0_funct(X,Y);  %1/(2*pi^2)*sin(x+y).^2;


%append ghost nodes around matrices

f1g = [f1(end,end),f1(end,:),f1(end,1);
       f1(:,end),f1, f1(:,1);
       f1(1,end), f1(1,:), f1(1,1)];
          
f2g = [f2(end,end),f2(end,:),f2(end,1);
       f2(:,end),f2, f2(:,1);
       f2(1,end), f2(1,:), f2(1,1)];   


[U_fd,T_fd] = AB2_method_2d(@(u,t) fd2(u,t,f1g,f2g,h),u0,dt,nsteps,iosteps);
         

end

%-----------------------------------------------------------------
%-- computes RHS of semi-discrete form of PDE --

function rhs = fd2(u,ts,f1,f2,h)

%[n,m] = size(u);


 ug = [u(end,end),u(end,:),u(end,1);
       u(:,end),u, u(:,1);
       u(1,end), u(1,:), u(1,1)];

           
rhs = zeros(size(ug)); %u is the ghost padded matrix

for j = 2:size(ug,2) - 1
    for i = 2:size(ug,1) - 1
        rhs(i,j) = -(f1(i+1,j)*ug(i+1,j) - f1(i-1,j)*ug(i-1,j))/(2*h) - ...
                    (f2(i,j+1)*ug(i,j+1) - f2(i,j-1)*ug(i,j-1))/(2*h);
    end
end

rhs = rhs(2:end-1,2:end-1); %make result back to normal size, no ghost padding


end
           
              
          

