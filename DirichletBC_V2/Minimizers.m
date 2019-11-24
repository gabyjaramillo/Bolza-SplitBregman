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

lmb = 0.01; %For constraint (appears in Gauss Seidel iteration)
h =0.01;


example = 'double'; % options are: 'double', 'double-half', 'triple'
potential = 'non-convex';  % can take values 'non-convex', 'convex'

% If picking potentia= non-convex, you have to pick value of g in section
% Matrices for Gauss-Seidel.

nmx= 2^8;                       % number of nodes
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

a = 4*abs(min(g))+0.1;

e = ones(nmx-2,1);
Upr = (lmb/dx^2)*spdiags(-e,1,nmx-2,nmx-2);

switch potential
    case 'non-convex'
        Lwr = (lmb/dx^2)*spdiags([-e 2*e],-1:0,nmx-2,nmx-2) +...
            (1/h+4*a)*speye(nmx-2); 
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

[u,error, count, Energy ,E_tot] = Split_Bregman_Combined(parameters, vals, dd, g, Lwr, Upr,u0,example,potential);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %   Plots and Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
plot(xx,u,'LineWidth',2)


M = count-2;
iter = 1:1:count-1;
figure(2)
semilogy(iter(end-M:end),Energy(end-M:end),'LineWidth',2)
set(gca,'FontSize',16)
set(0,'defaulttextInterpreter','latex')
xlabel('Iterations')
ylabel('$log_{10}(\bar{I}$)')

axes('Position',[0.4 0.4 .4 .4] )
box on
semilogy(iter(end-4000:end),Energy(end-4000:end),'LineWidth',2)
set(gca,'FontSize',14)


figure(3)
semilogy(iter(end-M:end),error(count-M:count),'LineWidth',2)
set(gca,'FontSize',16)
set(0,'defaulttextInterpreter','latex')
xlabel('Iterations')
ylabel('$\log_10$(Error)')
axes('Position',[0.55 0.55 .3 .3] )
box on
semilogy(iter(end-500:end),error(count-500:count),'LineWidth',2)
set(gca,'FontSize',14)


figure(4)
semilogy(iter(end-M:end),E_tot(end-M:end),'LineWidth',2)
set(gca,'FontSize',16)
set(0,'defaulttextInterpreter','latex')
xlabel('Iterations')
ylabel('$log_{10}(E_{grad}$)')

axes('Position',[0.4 0.4 .4 .4] )
box on
semilogy(iter(end-4000:end),E_tot(end-4000:end),'LineWidth',2)
set(gca,'FontSize',14)





% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %  Shanks Transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%S = (An+1 An-1 - An^2)/ (An+1 -2 An +An-1)
M=1000;
S = (Energy(end-3*M+1:end-2*M).*Energy(end-3*M-1:end-2-2*M) - Energy(end-3*M:end-1-2*M).^2)./...
(Energy(end-3*M+1:end-2*M) -2*Energy(end-3*M:end-1-2*M) + Energy(end-3*M-1:end-2-2*M) );


figure(5)
plot(S,'LineWidth',2)
set(gca,'FontSize',16)
set(0,'defaulttextInterpreter','latex')
title('Limiting Energy with Shanks')




%-----Extra ----

