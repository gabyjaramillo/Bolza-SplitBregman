%% Obstacle Problem (Chebishev points)
% The following code computes the solutions to the obstacle problem
%  J(u) = (1/2) int_Omega | \grad u|^2 dx 
% with dirichlet b.c. ( u \in H^1_0(Omega) ) and u >= W(x)
% Omega = [-1,1] and W(x) = (x^2-1)^2

%% Theory: from Tran and Osher
% Consider the functional 
% Jmu(u) = \int_Omega (1/2)|\grad u|^2 + mu (W-u)+ dx
% If u and u_mu minimize J(u) and Jmu(u), respectively. Then for mu >= -Wxx
% we have that u = u_mu

function uu = Obstacle(N,x,dx)
%=========================================================================
%Set up
% N =70;
% a=-1e-2;
% 
% b=2;
% dx =(b-a)/(N+1);
% x = (a:dx:b)';

mu =1e4;
lmbd=1;

Delta =  (-2*eye(N) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1))/(dx^2);
A = eye(N) - Delta/lmbd ;

W = 9-(x.^2-1).^2;
d = zeros(N,1);
b = d;
u = d;
%BC
u0 = W(1);
uend = W(end);
count = 1;
Tcount = 1e10;
tol = 1e-8;
err=1;

%=========================================================================
%=========================================================================
%Begin iteration

while (err>tol) && (count<Tcount)

%=========================================================================
%Gradient Descent: minimizing u
    uold =u;
    frhs = W(2:end-1) -d - b +(1/lmbd)*[u0/(dx^2); zeros(N-2,1);uend/(dx^2)];
    u = A\frhs;

%=========================================================================
% Shrink operator: minimizing d
% Shrink(z,c) = z-c if z>c , z if z<0, 0 otherwise

    z = W(2:end-1) - u - b;
    c = mu/lmbd;
    C = ones(N,1)*c;

    ind1 = find(z<0);
    ind2 = find(z>c);

    d = zeros(N,1);
    d(ind1) = z(ind1);
    d(ind2) = z(ind2) - C(ind2);

%=========================================================================
% b update

    b = b+u+d - W(2:end-1);

%=========================================================================
% Error
    unew = u;
    err = norm(unew-uold);
    error(count) = err;
    count = count +1;
    
end
u= [u0; u; uend];
uu = 9-u;

%  Waux = csapi(x,uu);
%   
%  Wrelax = fnval(Waux,x);
% 
% dW =( Wrelax(2:end)-Wrelax(1:end-1))/dx;
% disp(dW(8)+sqrt(32/27))

%  figure(10)
%  plot(x,Wrelax,'o-')
%  hold on
%  plot(x(2:end),dW)
%  plot(x,-ones(size(x))*sqrt(32/27))
%  plot(x,(x.^2-1).^2)
% hold off
%  legend('Relaxed W','line','(x^2-1)^2')
% 

% figure(2)
% plot(uu)
% hold on
% plot(9-W)
% hold off



% 
% figure(1)
% plot(x, u,'LineWidth',2)
% hold on
% plot(x,W)
% hold off

% figure(2)
% plot(x,uu,'o')
% hold on
% plot(x, Wrelax)
% hold off

