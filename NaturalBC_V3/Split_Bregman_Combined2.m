%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Split Bregman and Gradient Flow for Nonsmooth Nonconvex functionals
% Shankar Venkataramani & Gabriela Jaramillo

% E = int (u-u_k)^2/(2h) +W[u'] +V[u] dx  x in D

% Here W[u'] can be chosen to be the convex envelope of
% W1(d) = (d^2-1)^2 
% W2(d) = (d^2-1)^2 if d >= 0 and infty if d<0 
% W3(d) = (d^2-1)^2( (d-1)^2-1)^2 

% Also V[u] can be chosen to be
% V1[u] = (u-g)^2   
% V3[u] = (u^2 -g(x))^2 The function g(x) has to be >0 in this case.

% We use natural BC
% We will assume that the function lives on the nodes, the derivative
% lives on the intervals between the nodes.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u, error, count,Energy, U_min, E_min, count_min] = Split_Bregman_Combined2(par, vals, dd, g, Lwr, Upr,u0, example,potential)
global lmb a h

nmx = par(1);
dx  = par(2);
coef = par(3);

tol = 1.0e-12;     % tolerance

% Constants
nstep = 50000;  % maximum number of iterations 
nGS = 20;       % number of Gauss-Seidel iterations per update of u,v and b.
kk = 5;         % number of Gradient Flow iterations
error = zeros(1,nstep);
Energy = error;
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
U_min = zeros(size(u0));
E_min = 100;


while((err > tol) && (count <= nstep))
    
    uold = u;
    for ii=1:kk  %Gradient Flow iteration 
    
% Gauss-Siedel update for u  
%---------------------------
    w = v-b;
    wext = [-w(1); w ; -w(end)]; %This accounts for natural BC
 
    for jj=1:nGS
        frcng =  - Upr*u - uold./h + lmb*(wext(2:end) - wext(1:end-1))./dx ...       
                 +coef*(-4*(a+g).*uold + 4*uold.^3) + (1-coef)*(-2*g);
             
       
        u = Lwr\frcng;
    end

 
% Shrink Operator update for v
%-------------------------------  
    vtmp = (u(2:end,1)-u(1:end-1,1) )/dx +b;
    v = interp1(z_shrink, d_shrink, vtmp);

% Update for b
%--------------
    b = b-v+(u(2:nmx)-u(1:nmx-1))/dx;   
    
    end
    
    % Calculating Energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ux = ( u(2:end)-u(1:end-1) )./ dx;
        
switch example
    case 'double'
        E1 =(ux.^2 - 1).^2;
        ind = find(abs(ux)<=1);
        E1(ind) = 0;
        
    case 'double-half'
        E1 = (ux.^2 - 1).^2;
        
    case 'triple'
        E1 = (ux.^2 - 1).^2.*((ux-2).^2-1).^2;
        ind = find(-1<= ux<=3);
        E1(ind) = 0;
        
end

switch potential
    case 'convex'
        E2 = (u - g).^2;
    case 'non-convex'
        E2 = (u.^2-g).^2;
        
end

E1 = sum((E1(1:end-1) + E1(2:end))*(dx/2));
E2 = sum((E2(1:end-1) + E2(2:end))*(dx/2));

EE = E1 +E2;
Energy(count) = EE;

if EE < E_min
    U_min = u;
    E_min = EE;
    count_min = count;
end

    
    
% Checking error   
%-----------------    
    err = norm(msk.*(u-uold))*sqrt(dx);
    error(count) = err;
    count = count+1;
    
end

end


