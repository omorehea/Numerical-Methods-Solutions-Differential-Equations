function [x,w,D] = gen_GLL_points_and_D_matrix(N,a,b)

% This function computes the Gauss-Legendre-Lobatto
% differentiation matrix in the interval [a,b]
%
% Input N   -> polynomial order
%       a,b -> interval [a,b]
%
% Output x  -> vector of (N+1) GLL quadrature points  in [a,b].
%        w  -> vector of (N+1) GLL quadrature weights in [a,b].
%        D  -> (N+1) x (N+1) diff matrix such that df/dx = D*f in [a,b]
%



N1=N+1;

% Chebyshev Gauss Lobatto nodes
xc=cos(pi*(0:N)/N)';

% Uniform nodes
xu=linspace(-1,1,N1)';

% Make a close first guess to reduce iterations
if N<3
    x=xc;
else
    x=xc+sin(pi*xu)./(4*N);
end

% Legendre Vandermonde Matrix
P=zeros(N1,N1);

% Compute P(N) using the recursion relation
% Compute its first and second derivatives and 
% update x using the Newton-Raphson method.

xold=2;
while max(abs(x-xold))>eps

    xold=x;
        
    P(:,1)=1;    P(:,2)=x;
    
    for k=2:N
        P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
    end
     
    x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
end

X=repmat(x,1,N1);
Xdiff=X-X'+eye(N1);

L=repmat(P(:,N1),1,N1);
L(1:(N1+1):N1*N1)=1;
D=(L./(Xdiff.*L'));
D(1:(N1+1):N1*N1)=0;
D(1)=(N1*N)/4;
D(N1*N1)=-(N1*N)/4;

w=2./(N*N1*P(:,N1).^2);

x    = flipud(x);
w    = flipud(w);
D    = -D;

D=2/(b-a)*D;
w=(b-a)/2*w;
x=(b-a)/2*x+(b+a)/2;

