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
clear; clc; close all
global lmb a h

%a = 4;   %For convex splitting, a>2
%lmb = 0.01; %For constraint (appears in Gauss Seidel iteration)
%h =0.01;

example = 'double-half'; % options are: 'double', 'double-half', 'triple'
potential = 'convex';  % can take values 'non-convex', 'convex'

alpha = 0; beta =1;  %end points of x interval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Shrink Operator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=150; % Number of grid points for Obstacle problem

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
    
nmx = [2^5 2^6 2^7 2^8 2^9 2^10];
T = zeros(size(nmx));
Energy = T;
for ii = 1: length(nmx)
h = max(1./nmx(ii),0.01);
 lmb = h;

% Possible values of g(x)
%-------------------------------------------
           % number of nodes
dx = (beta- alpha)/(nmx(ii)-1);     % grid spacing
xx = (alpha:dx:beta)'; 


% g = sin(2*pi*xx)/4;
 g = zeros(size(xx));
% g = -1/2*xx;
% g = sin(2*pi*xx)/6+exp(xx)/2;
% g = exp(xx);
% g = sin(4*pi*xx)/12+exp(xx)/2;
% g = -(3/128*16)*(xx).^5 - (xx).^3/12;
% g = -(128/3)*(abs(xx)-0.5).^5 - (abs(xx)-0.5).^3/3;
% g = 1.5*xx;

a = 4*abs(min(g))+1;
e = ones(nmx(ii)-2,1);
Upr = (lmb/dx^2)*spdiags(-e,1,nmx(ii)-2,nmx(ii)-2);


switch potential
    case 'non-convex'
        Lwr = (lmb/dx^2)*spdiags([-e 2*e],-1:0,nmx(ii)-2,nmx(ii)-2) +...
            (1/h+4*a)*speye(nmx(ii)-2); 
        coef = 1;
        
    case 'convex'
        Lwr = (lmb/dx^2)*spdiags([-e 2*e],-1:0,nmx(ii)-2,nmx(ii)-2) + ...
            (2+ 1/h)*speye(nmx(ii)-2); 
        coef = 0;
      
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Call on Modified Split Bregman
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u0 = 0.1*ones(size(xx));  %include desired boundary conditions
u_L=0;
u_R=0.5;

parameters = [nmx(ii),dx,u_L,u_R, coef];
tic
u = Split_Bregman_Combined(parameters, vals, dd, g, Lwr, Upr,u0,example,potential);
T(ii) = toc;


figure(3)
 plot(xx,u, 'LineWidth',2)
hold on
end
hold off

