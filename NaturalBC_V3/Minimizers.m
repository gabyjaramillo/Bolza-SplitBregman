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
% the function g(x) can be chose in section "Matrices for Gauss-Seidel"

% We assume natural BC ( u_x = b+d) at end points

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

nmx= 2^8;                       % number of nodes
alpha = -1; beta =1;             % end points
dx = (beta-alpha)/(nmx-1);      % grid spacing
xx = (alpha:dx:beta)'; 

u0 = 0.1*ones(size(xx));    % initial guess and BC


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Shrink Operator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=200; % Number of grid points for Obstacle problem

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
 g = sin(2*pi*xx)/4 + 1/2;
% g = sin(2*pi*xx)/4;
% g = ones(size(xx));
% g = sin(2*pi*xx)/6+exp(xx)/2;
% g = exp(xx);
% g = sin(4*pi*xx)/12+exp(xx)/2;
% g = -(128/3)*(abs(xx)-0.5).^5 - (abs(xx)-0.5).^3/3;
% g = 1.5*xx;

a = 4*abs(min(g))+0.1; %for convex splitting (u^2-g)^2 = 2au^2 +(u^4- 2(g+a)u^2 + g^2)
e = ones(nmx,1);
Upr = (lmb/dx^2)*spdiags(e,1,nmx,nmx);
Upr(1,2) = (lmb/dx^2)*(2);


switch potential
    case 'non-convex'
        Lwr = (lmb/dx^2)*spdiags([e -2*e],-1:0,nmx,nmx)...
            - (1/h+4*a)*speye(nmx);
        Lwr(end,end-1) = (lmb/dx^2)*(2);
        coef = 1;
        
    case 'convex'
        Lwr = (lmb/dx^2)*spdiags([e -2*e],-1:0,nmx,nmx) - (1/h +2)*speye(nmx);
        Lwr(end,end-1) = (lmb/dx^2)*(2);
        coef = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Call on Modified Split Bregman
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


parameters = [nmx,dx, coef];

[u,error, count,Energy, U_min, E_min, count_min] = Split_Bregman_Combined2(parameters, vals, dd, g, Lwr, Upr,u0,example,potential);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %   Plots and Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Ux = (U_min(2:end) - U_min(1:end-1))./dx;

%Plot u and ux
figure(1) 
[hAx,hLine1,hLine2] = plotyy(xx,U_min,xx(1:end-1),Ux);
set(groot, 'defaultLegendInterpreter','latex');
lg = legend('$u(x)$','$u_x(x)$');
hLine2.LineStyle = '--';
hLine1.LineWidth = 2;
hLine2.LineWidth = 2;
hAx(1).FontSize = 16;
hAx(2).FontSize = 16;
hAx(1).YColor='k';
hAx(2).YColor ='k';
hAx(1).YGrid ='on';
hAx(1).XGrid ='on';
%hAx(1).XTick =[-1 -0.66 0 1]; % This needs to change for
%different graphs
xlabel('$x$')
ylabel(hAx(1),'$u(x)$') % left y-axis 
ylabel(hAx(2),'$u_x(x)$') % right y-axis
lg.FontSize =16;



%%% For non-convex
%Plot ux and derivative of sqrt(g)

d_g = (g(2:end)-g(1:end-1))/dx;
dd_sqrtg = (1/2)*d_g./(sqrt(g(2:end)));

figure(2)
plot(xx(1:end-1), Ux, 'LineWidth',2)
set( gca,'FontSize',16)
set(groot, 'defaultLegendInterpreter','latex');
hold on
plot(xx(1:end-1), dd_sqrtg ,'--', 'LineWidth',2)
hold off
lg =legend('$u_x$','$(\sqrt{g})_x$');
lg.FontSize = 16;
title('Nonconvex')


%Plot u and sqrt(g)
figure(3)
plot(xx, sqrt(g) , 'LineWidth',2)
set( gca,'FontSize',16)
set(groot, 'defaultLegendInterpreter','latex');
hold on
plot(xx(1:4:end), U_min(1:4:end),'*', 'LineWidth',2)
hold off
lg =legend('$\sqrt{g(x)}$','u(x)');
lg.FontSize = 16;
xlabel('$x$')
%title('Nonconvex')


%%% For convex potential
%%Plot of ux and derivative of g

figure(4)
plot(xx(1:end-1), Ux, 'LineWidth',2)
set( gca,'FontSize',16)
set(groot, 'defaultLegendInterpreter','latex');
hold on
plot(xx(1:end-1), d_g ,'--', 'LineWidth',2)
hold off
lg =legend('$u_x$','$g_x$');
lg.FontSize = 16;
title('Convex')

%Plot u and g
figure(5)
plot(xx, g , 'LineWidth',2)
set( gca,'FontSize',16)
set(groot, 'defaultLegendInterpreter','latex');
hold on
plot(xx(1:8:end), U_min(1:8:end),'*-', 'LineWidth',1)
hold off
lg =legend('$g(x)$','$u(x)$');
lg.FontSize = 16;
xlabel('$x$')
%title('Convex')



%Plot difference between limiting solution, u, and minimizer, U_min

figure(6)
plot(xx, abs(U_min -u),'LineWidth',2)
set( gca,'FontSize',16)
set(groot, 'defaultLegendInterpreter','latex');
xlabel('$x$')
ylabel('$|U_{min} -u|$')


M = count-2;
iter = 1:1:count-1;

%Plots energy at each iteration of Gradient Flow and specifies the
%iteration with the lowest energy
figure(7)
clf
semilogy(iter,Energy(1:count-1),'LineWidth',2)
set(gca,'FontSize',16)
set(0,'defaulttextInterpreter','latex')
hold on 
semilogy(count_min,E_min,'*','LineWidth',2)
xlabel('Iterations')
ylabel('log scale')
title('Energy $\bar{I}$')
%  axes('Position',[0.4 0.4 .4 .4] )
%  box on
%  semilogy(iter(1400:2600),Energy(1400:2600),'LineWidth',2)
%  set(gca,'FontSize',14)
%  hold on
%  semilogy(count_min,E_min,'*','LineWidth',2)
%  semilogy(iter(1400:2600),ones(2600-1399,1)*1.02408,'LineWidth',2)
%  hold off
%  axis([ 1400 2600 1.0237788 1.02379])


%Plots the relative error, |u_n - u_n+1|, at each iteration of Gradient
%Flow
figure(8)
semilogy(iter,error(1:count-1),'LineWidth',2)
set(gca,'FontSize',16)
set(0,'defaulttextInterpreter','latex')
xlabel('Iterations')
ylabel('Log scale')
%axes('Position',[0.55 0.55 .3 .3] )
% box on
% semilogy(iter(end-500:end),error(count-500:count),'LineWidth',2)
% set(gca,'FontSize',14)
