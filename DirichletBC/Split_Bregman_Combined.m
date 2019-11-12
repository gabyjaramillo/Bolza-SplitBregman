%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Split Bregman and Gradient Flow for Nonsmooth Nonconvex functionals
% Shankar Venkataramani & Gabriela Jaramillo

% E = int (u-u_k)^2/(2h) +W[u'] +V[u] dx  x in D

% Here W[u'] can be chosen to be the convex envelope of
% W1(d) = (d^2-1)^2 
% W2(d) = (d^2-1)^2 if d >= 0 and infty if d<0 
% W3(d) = (d^2-1)^2( (d-1)^2-1)^2 

% Also V[u] can be chosen to be
% V1[u] = u^2   
% V2[u] = (u^2-1)^2
% V3[u] = (u^2 -g(x))^2 


% We will assume that the function lives on the nodes, the derivative
% lives on the intervals between the nodes.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u, error, count] = Split_Bregman_Combined(par, vals, dd, g, Lwr, Upr,u0)
%u = Split_Bregman_Combined(nmx,dx, vals,dd, g, Lwr, Upr,coef, u0)
global lmb a h

nmx = par(1);
dx  = par(2);
u_L = par(3);
u_R = par(4);
coef = par(5);

tol = 1.0e-12;     % tolerance

% Constants
nstep = 30000;  % maximum number of iterations 
nGS = 10;       % number of Gauss-Seidel iterations per update of u,v and b.
kk = 2;         % number of Gradient Flow iterations
error = zeros(1,nstep);
err = 1.0;


% weights to do trapezoid rule for computing the error
%--------------------------------------------------------
msk = ones(nmx,1); 
msk(1) = 1/sqrt(2.);
msk(end) = 1/sqrt(2.);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Shrink Operator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  We compute the subgradients at the "breaks" to identify the inverse of
%  the shrink operator


slopes = (vals(2:end)-vals(1:end-1))./(dd(2:end)-dd(1:end-1));
curvature  = (slopes(2:end)-slopes(1:end-1))./(dd(2)-dd(1));

ind = find(curvature >1e-6); % find discontinuities in D(W[d]) 
auxind= [1; ind+1];
dd_aux= dd(auxind)'; %because curvature is offset by one in the index

vals2= vals(auxind)';
slopes2 = (vals2(2:end)-vals2(1:end-1))./(dd_aux(2:end)-dd_aux(1:end-1));

subgrad = [-1000 slopes2 1000;-1000 slopes2 1000 ];
subgrad = subgrad(:)';

d_shrink = [dd_aux;dd_aux];
d_shrink =[-100 d_shrink(:)' 100];

z_shrink = d_shrink + subgrad/lmb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Iteration Modified Split Bregman
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initializing variables
b = zeros(nmx-1,1);
v = b;
u = u0;
count = 1; % count number of iterations


while((err > tol) && (count <= nstep))
    
    uold = u;
    for ii=1:kk  %Gradient Flow iteration 
    
% Gauss-Siedel update for u  
%---------------------------
    wext = v-b;
    u(1) = u_L; u(end) = u_R;  % Boundary Conditions
    
    for jj=1:nGS
        frcng = -Upr*u(2:end-1)+(lmb/dx^2)*[u(1);zeros(nmx-4,1);u(end)]...
              + lmb*(wext(1:nmx-2)-wext(2:nmx-1))/dx + (1/h)*uold(2:end-1)...
              + coef*(4*(1+a)*g(2:end-1).*uold(2:end-1) - 4*uold(2:end-1).^3);
          
        u(2:end-1) = Lwr\frcng;
    end

 
% Shrink Operator update for v
%-------------------------------  
    vtmp = (u(2:end,1)-u(1:end-1,1) )/dx +b;
    v = interp1(z_shrink, d_shrink, vtmp);

% Update for b
%--------------
    b = b-v+(u(2:nmx)-u(1:nmx-1))/dx;   
    
    end
    
% Checking error   
%-----------------    
    err = norm(msk.*(u-uold))*sqrt(dx);
    error(count) = err;
    count = count+1;
    
end
disp(err)
disp(count)
end


