function [x,D1,D2] = get_GCL_points_and_D_matrices(a,b,N)

% This function returns N+1 Gauss Chebyshev Lobatto (GCL)
% point in [a,b] and the corresponding differentiation 
% matrices D1 (first-order) and D2 (second-order).
%
% To compute the first and second derivative of a function f(x) at 
% the GCL nodes x it is sufficient to compute
% D1*f(x) and D2*f(x)
%
% Example of application:
% 
% [x,D1,D2] = get_GCL_points_and_D_matrix(0,2*pi,100);
%
% figure(1) 
% plot(x,sin(2*x))
% hold
% plot(x,D1*sin(2*x))     % first  derivative of sin(2*x) =  2*cos(2*x)
% plot(x,D2*sin(2*x))     % second derivative of sin(2*x) = -4*sin(2*x)



eta=fliplr(cos(pi/N*(0:N)))'; % GCL points in [-1,1]

x = (b-a)/2*eta + (b+a)/2;


% Here we build the first-order differentiation matrix
d=ones(1,N+1);
d(1)=2; 
d(end)=2;

D1=zeros(N+1,N+1);
D2=zeros(N+1,N+1);
D1 = D1 + diag(-eta./(2*(1-eta.^2)));

for ii=0:N
    for jj=0:N
        if ii~=jj
             D1(ii+1,jj+1)= d(ii+1)/d(jj+1) * (-1)^(ii+jj)/(eta(ii+1)- eta(jj+1));
        end
        
    end
end

D1(1,1)    = -(2*N^2+1)/6; 
D1(end,end)= -D1(1,1);




% Here we build the second-order differentiation matrix
eta=cos(pi/N*(0:N))';

for k=1:N-1
    for l=0:N
        if k~=l
             D2(k+1,l+1)= (-1)^(k+l) * ...
                           (  eta(k+1)^2 + eta(k+1)*eta(l+1) -2) / ... 
                           (  (1 - eta(k+1)^2) * (eta(k+1) - eta(l+1))^2  )/d(l+1);
        end
        
    end
end

for ii=1:N-1
            D2(ii+1,ii+1)= - ((N^2-1)*(1-eta(ii+1)^2)+3) / (3*(1-eta(ii+1)^2)^2); 
end


for jj=1:N
            D2(1,jj+1)= 2*(-1)^(jj)/(3*d(jj+1))*  ((2*N^2+1)*(1-eta(jj+1))-6)/(1-eta(jj+1))^2; 
           
end

for jj=0:N-1             
            D2(N+1,jj+1)= 2*(-1)^(jj)/(3*d(jj+1))*((2*N^2+1)*(1+eta(jj+1))-6)/(1+eta(jj+1))^2;   
end

D2(1,1)= (N^4-1)/15;
D2(end,end)= (N^4-1)/15;


D2=rot90(D2,2);

D1  = 2/(b-a)*D1;
D2 = (2/(b-a))^2*D2;










end

