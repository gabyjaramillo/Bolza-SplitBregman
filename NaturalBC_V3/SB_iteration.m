% Numerically find minimizers for Nonsmooth potentials using Modified 
% Split Bregman 

% Gabriela Jaramillo & Shankar Venkataramani


% E = int (u-u_k)^2/(2h) +W[u'] +V[u] dx  x in D

% Here W[u'] can be chosen to be the convex envelope of
% W1(d) = (d^2-1)^2  "double"
% W2(d) = (d^2-1)^2 if d >= 0 and infty if d<0     "double-half"
% W3(d) = (d^2-1)^2( (d-1)^2-1)^2   "triple"

% Also V[u] can be chosen to be
% V1[u] = (u- g(x))^2   "convex potential"
% V2[u] = (u^2 -g(x))^2  "non-convex potential"

% We assume natural BC

% We compute the convex envelop of W(d) as an obstacle problem

% We will assume that the function lives on the nodes, the derivative
% lives on the intervals between the nodes.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all
global lmb a h


example = 'double'; % options are: 'double', 'double-half', 'triple'
potential = 'non-convex';  % can take values 'non-convex', 'convex'

alpha = -1; beta =1;  %end points of x interval
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


 g = sin(2*pi*xx)/4 +1/2;
 %g = sin(2*pi*xx)/4 +1/2;
% g = ones(size(xx));
% g = -1/2*xx;
 %g = sin(2*pi*xx)/6+exp(xx)/2;
% g = exp(xx);
% g = sin(4*pi*xx)/12+exp(xx)/2;
% g = -(3/128*16)*(xx).^5 - (xx).^3/12;
% g = -(128/3)*(abs(xx)-0.5).^5 - (abs(xx)-0.5).^3/3;
% g = 1.5*xx;

a = 4*abs(min(g))+0.1;

e = ones(nmx(ii),1);
Upr = (lmb/dx^2)*spdiags(e,1,nmx(ii),nmx(ii));
Upr(1,2) = (lmb/dx^2)*(2);


switch potential
    case 'non-convex'
        Lwr = (lmb/dx^2)*spdiags([e -2*e],-1:0,nmx(ii),nmx(ii))...
            - (1/h+ 4*a)*speye(nmx(ii));
        Lwr(end,end-1) = (lmb/dx^2)*(2);
        coef = 1;
        
    case 'convex'
        Lwr = (lmb/dx^2)*spdiags([e -2*e],-1:0,nmx(ii),nmx(ii)) - (1/h +2)*speye(nmx(ii));
        Lwr(end,end-1) = (lmb/dx^2)*(2);
        coef = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Call on Modified Split Bregman
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u0 = 0.1*ones(size(xx));  %include desired boundary conditions


parameters = [nmx(ii),dx coef];
tic
[u, error, count,Energy, U_min, E_min, count_min] = Split_Bregman_Combined2(parameters, vals, dd, g, Lwr, Upr,u0, example, potential);
T(ii) = toc;

E(ii) = E_min;
Error(ii) = error(count-1);
count_iter(ii) = count;
E2(ii) = Energy(count-1);

U_x = (U_min(2:end) - U_min(1:end-1))./dx;
u_x = (u(2:end) -u(1:end-1))./dx;

% figure(1)
% plot(xx(1:end-1),U_x,'LineWidth',2)
% hold on
% plot(xx(1:end-1),u_x,'--','LineWidth',1)

figure(1)
plot(xx, abs(u-U_min))
hold on

end
hold off

