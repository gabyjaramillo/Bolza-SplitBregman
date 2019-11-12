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

function u = Obstacle(string,N,x,dx)
%=========================================================================


mu =1e4;
lmbd=1;

Delta =  (-2*eye(N) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1))/(dx^2);
A = eye(N) - Delta/lmbd ;

W = string(x);
d = zeros(N,1);
b = d; u = d;

u0 = W(1); uend = W(end); %BC
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
    count = count +1;
    
end
u= [u0; u; uend];


