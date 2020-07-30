% Numerically find minimizers for Nonsmooth potentials using Modified 
% Split Bregman 

% Gabriela Jaramillo & Shankar Venkataramani


% E = int (u-u_k)^2/(2h) +W[u'] +V[u] dx  x in D

% Here W[u'] can be chosen to be the convex envelope of
% W1(d) = (d^2-1)^2  "double"
% W2(d) = (d^2-1)^2 if d >= 0 and infty if d<0     "double-half"
% W3(d) = (d^2-1)^2( (d-1)^2-1)^2   "triple"

% Also V[u] can be chosen to be
% V1[u] = (u-g(x))^2  "convex potential"
% V2[u] = (u^2 -g(x))^2 "non-convex potential"

% We compute the convex envelop of W(d) as an obstacle problem

% We will assume that the function lives on the nodes, the derivative
% lives on the intervals between the nodes.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all
global lmb a h

potential = 'convex';  % can take values 'non-convex', 'convex'
example = 'double-half'; % options are: 'double', 'double-half', 'triple'

%if 'random' example is selected then:
from_file = 'no'; %options are: 'yes' if random functions is from file, 'no' if random function is defined here

alpha =0; beta =1;  %end points of x interval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Convex Hull using Beneath and Beyond Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=400; % Number of grid points for convex hull algorithm
 
switch example
    case 'double'
        a0 = -2;
        b0 = 2;
        deltad =(b0-a0)/(N+1);
        d = (a0:deltad:b0)';
        well = (d.^2-1).^2;
        sm = -100;
        sp = 100;
        
    case 'double-half'
        
        a0 = 0;
        b0 = 2;
        deltad =(b0-a0)/(N+1);
        d = (a0:deltad:b0)';  
        well = (d.^2-1).^2;
        sm = -100;
        sp = 100;
        
    case 'triple'
       
        a0 = -2;
        b0 = 4;
        deltad =(b0-a0)/(N+1);
        d = (a0:deltad:b0)'; 
        well = (d.^2-1).^2.*( (d-2).^2 -1).^2;
        sm = -10^4;
        sp = 10^4;
        
    case 'random'
        aux_2 = strcmp(from_file,'yes');        
        
        if aux_2 ==1
            fileID = fopen('random_potential.txt','r');
            formatSpec = '%f';
            A = fscanf(fileID,formatSpec);
            [aux,~] = size(A);
            B = reshape(A,[aux/2,2]);
            d = B(:,1)';
            well = B(:,2)';
            a0 = d(1);
            b0 = d(end);
            sm = -5/4;
            sp = 5/4;
            
        else
            a0 = -2;
            b0 = 2;
            sm = -5*rand/2;
            sp = 5*rand/2;
            d = linspace(a0,b0,11);
            well = rand(1,11);
            
        end
end

%--bb method
[CX,CY] = bb(d,well,sm,sp);
aux = strcmp(example,'random');

if aux ==1;
    
    dd = [a0-1 CX b0+1]';
    vals = [well(1)-sm, CY ,well(end)+sp]';
 
    d = [a0-1 d b0+1]';
    well = [well(1)-sm well well(end)+sp]';
    

else
    dd = CX';
    vals = CY';
   
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Matrices for Guass Seidel with Dirichlet BC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
nmx = [2^5 2^6 2^7 2^8 2^9 2^10];
T = zeros(size(nmx));
E = T;
E2 = T;
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
u_R=1/2;

parameters = [nmx(ii),dx,u_L,u_R, coef];
tic
[u, count, Energy, U_min, E_min, count_min,error_H_count] = Split_Bregman_Combined(parameters, vals, dd, g, Lwr, Upr,u0,example,potential);
T(ii) = toc;
E(ii)= E_min;
E2(ii) = Energy(count-1);


end




