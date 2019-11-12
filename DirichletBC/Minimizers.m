% Numerically find minimizers for Nonsmooth potentials using Modified 
% Split Bregman 

% Gabriela Jaramillo & Shankar Venkataramani


% E = int (u-u_k)^2/(2h) +W[u'] +V[u] dx  x in D

% Here W[u'] can be chosen to be the convex envelope of
% W1(d) = (d^2-1)^2  "double"
% W2(d) = (d^2-1)^2 if d >= 0 and infty if d<0     "double-half"
% W3(d) = (d^2-1)^2( (d-1)^2-1)^2   "triple"

% Also V[u] can be chosen to be
% V1[u] = u^2   
% V2[u] = (u^2-1)^2
% V3[u] = (u^2 -g(x))^2 

% We always assume homogeneous Dirichlet BC

% We compute the convex envelop of W(d) as an obstacle problem

% We will assume that the function lives on the nodes, the derivative
% lives on the intervals between the nodes.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc;
global lmb a h

a = 3/2;   %For convex splitting, a>2
lmb = 2; %For constraint (appears in Gauss Seidel iteration)
h =0.01;


example = 'triple'; % options are: 'double', 'double-half', 'triple'
potential = 'non-convex';  % can take values 'non-convex', 'convex'
% If picking potentia= non-convex, you have to pick value of g in section
% Matrices for Gauss-Seidel.

nmx= 2^6;                       % number of nodes
alpha = -1; beta =1;             % end points
dx = (beta-alpha)/(nmx-1);      % grid spacing
xx = (alpha:dx:beta)'; 

u0 = 0.1*ones(size(xx));  %initial guess and BC
u_L=0; u_R=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Obstacle Problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=400; % Number of grid points for Obstacle problem

switch example
    case 'double'
        well = @(x) 9- (x.^2-1).^2;       
        a0 = -2;
        b0 = 2;
        deltad =(b0-a0)/(N+1);
        dd = (a0:deltad:b0)';
        offset = 9;
        
    case 'double-half'
        
        well = @(x) 9- (x.^2-1).^2;      
        a0 = 0;
        b0 = 2;
        deltad =(b0-a0)/(N+1);
        dd = (a0:deltad:b0)';        
        offset =9;
        
    case 'triple'
        
        well = @(x) 2025-(x.^2-1).^2.*( (x-2).^2 -1).^2;      
        a0 = -2;
        b0 = 4;
        deltad =(b0-a0)/(N+1);
        dd = (a0:deltad:b0)';       
        offset =2025;

end

vals = offset - Obstacle(well,N,dd,deltad);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Matrices for Guass Seidel with Dirichlet BC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Possible values of g(x)
%-------------------------------------------
% g = sin(2*pi*xx)/4;
 g = ones(size(xx));
% g = -1/2*xx;
% g = sin(2*pi*xx)/6+exp(xx)/2;
% g = exp(xx);
% g = sin(4*pi*xx)/12+exp(xx)/2;
% g = -(3/128*16)*(xx).^5 - (xx).^3/12;
% g = -(128/3)*(abs(xx)-0.5).^5 - (abs(xx)-0.5).^3/3;
% g = 1.5*xx;

e = ones(nmx-2,1);
Upr = (lmb/dx^2)*spdiags(-e,1,nmx-2,nmx-2);

switch potential
    case 'non-convex'
        Lwr = (lmb/dx^2)*spdiags([-e 2*e],-1:0,nmx-2,nmx-2) +...
            4*a.*spdiags(g(2:end-1),0,nmx-2,nmx-2) + 1/h*speye(nmx-2); 
        coef = 1;
    case 'convex'
        Lwr = (lmb/dx^2)*spdiags([-e 2*e],-1:0,nmx-2,nmx-2) + ...
            (2+ 1/h)*speye(nmx-2); 
        coef = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Call on Modified Split Bregman
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


parameters = [nmx,dx,u_L,u_R, coef];

[u,error, count] = Split_Bregman_Combined(parameters, vals, dd, g, Lwr, Upr,u0);


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ux = ( u(2:end)-u(1:end-1) )./ dx;
        
switch example
    case 'double'
        E1 = (ux.^2 - 1).^2;
    case 'double-half'
        E1 = (ux.^2 - 1).^2;
    case 'triple'
        E1 = (ux.^2 - 1).^2.*((ux-2).^2-1).^2;
end

switch potential
    case 'convex'
        E2 = u.^2;
    case 'non-convex'
        E2 = (u.^2-g).^2;
end

E1 = sum((E1(1:end-1) + E1(2:end))*(dx/2));
E2 = sum((E2(1:end-1) + E2(2:end))*(dx/2));

energy = E1 +E2;
fprintf('Energy is %f\n',energy)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %   Plots and Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(10)
plot(xx,u,'LineWidth',2)

