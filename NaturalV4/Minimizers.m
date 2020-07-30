% Numerically find minimizers for Nonsmooth potentials using Modified 
% Split Bregman 

% Gabriela Jaramillo & Shankar Venkataramani


% E = int (u-u_k)^2/(2h) +W[u'] +V[u] dx  x in D

% Here W[u'] can be chosen to be the convex envelope of
% W1(d) = (d^2-1)^2  "double"
% W2(d) = (d^2-1)^2 if d >= 0 and infty if d<0     "double-half"
% W3(d) = (d^2-1)^2( (d-1)^2-1)^2   "triple"
% W4(d) = random piecewise function

% Also V[u] can be chosen to be
% V1[u] = (u- g(x))^2   "convex potential"
% V2[u] = (u^2 -g(x))^2  "non-convex potential"
% the function g(x) can be chose in section "Matrices for Gauss-Seidel"

% We assume natural BC ( u_x = b+d) at end points

% We compute the convex envelop of W(d) using the 'Beneath and Bellow'
% algorithm from Y. Lucet

% We will assume that the function lives on the nodes, the derivative
% lives on the intervals between the nodes.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc;
global lmb a h

lmb = 0.01; %For constraint (appears in Gauss Seidel iteration)
h =0.01;

potential = 'non-convex';  % can take values 'non-convex', 'convex'
example = 'random'; % options are: 'double', 'double-half', 'triple', 'random'
from_file = 'yes'; %options are: 'yes' if random functions is from file, 'no' if random function is defined here

nmx= 2^8;                       % number of nodes
alpha = -1; beta =1;             % end points
dx = (beta-alpha)/(nmx-1);      % grid spacing
xx = (alpha:dx:beta)'; 

u0 = 0.5*ones(size(xx));    % initial guess and BC


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       ux Potential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=200; % Number of grid points for convex hull algorithm
 
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


%-- Plot of convex hull and potential well

  figure(1)
  plot(dd,vals,'-o','LineWidth',2)
  set( gca,'FontSize',16)
  hold on
  plot(d,well,'LineWidth',2)
  hold off
  xg =xlabel({'$d$'},'Interpreter','latex');
  xg.FontSize =16;
  lg =legend({'$\overline{W}(d)$','$W(d)$'},'Interpreter','latex');
  lg.FontSize = 16;
 axis([dd(1),dd(end),-0.05,3])
   xt={-3 -2 -1.6 0.4 1.6 3} ; 
set(gca,'xtick',[-3 -2 -1.6 0.4 1.6 3]); 
set(gca,'xticklabel',xt);


 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %             Matrices for Guass Seidel with Dirichlet BC
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
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
Upr(1,2) = (lmb/dx^2)*(2); %Neumann BC

 
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

[u,count,Energy, U_min, E_min, count_min, error_H_count] = Split_Bregman_Combined2(parameters, vals, dd, g, Lwr, Upr,u0, example,potential);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    %   Plots and Figures
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%-- Plot of limiting solution, u, and U_min


Ux = (U_min(2:end) - U_min(1:end-1))./dx;
curve = ( Ux(2:end) - Ux(1:end-1))./dx;
curve_ind = find( (curve< -10) | ( curve>10));
x_tick = xx(curve_ind);
y_tick = Ux(curve_ind);


%Plot u and ux
figure(2) 
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
hAx(2).YGrid ='on';
hAx(2).XGrid ='on';
hAx(2).YTick = [-1 -0.5 0 0.5 1];
%hAx(1).XTick =[-1 -0.66 0 1]; % This needs to change for
%different graphs
xlabel('$x$')
ylabel(hAx(1),'$u(x)$') % left y-axis 
ylabel(hAx(2),'$u_x(x)$') % right y-axis
lg.FontSize =16;

% %%% For non-convex

% %Plot ux and derivative of sqrt(g)
% 
% d_g = (g(2:end)-g(1:end-1))/dx;
% dd_sqrtg = (1/2)*d_g./(sqrt(g(2:end)));
% 
% % figure(3)
% plot(xx(1:end-1), Ux, 'LineWidth',2)
% set( gca,'FontSize',16)
% set(groot, 'defaultLegendInterpreter','latex');
% hold on
% plot(xx(1:end-1), dd_sqrtg ,'--', 'LineWidth',2)
% hold off
% lg =legend('$u_x$','$(\sqrt{g})_x$');
% lg.FontSize = 16;
% title('Nonconvex')


% %Plot u and sqrt(g)

figure(4)
plot(xx, sqrt(g) , 'LineWidth',2)
set( gca,'FontSize',16)
%set(groot, 'defaultLegendInterpreter','latex');
hold on
plot(xx(1:4:end), U_min(1:4:end),'*-')
hold off
lg =legend({'$\sqrt{g(x)}$','$u(x)$'},'Interpreter','latex');
lg.FontSize = 16;
xlabel({'$x$'},'Interpreter','latex')
%title('Nonconvex')


% %%% For convex potential

% %%Plot of ux and derivative of g
% 
% figure(5)
% plot(xx(1:end-1), Ux, 'LineWidth',2)
% set( gca,'FontSize',16)
% set(groot, 'defaultLegendInterpreter','latex');
% hold on
% plot(xx(1:end-1), d_g ,'--', 'LineWidth',2)
% hold off
% lg =legend('$u_x$','$g_x$');
% lg.FontSize = 16;
% title('Convex')
% 
%Plot u and g

% figure(6)
% plot(xx, g , 'LineWidth',2)
% set( gca,'FontSize',16)
% set(groot, 'defaultLegendInterpreter','latex');
% hold on
% plot(xx(1:8:end), U_min(1:8:end),'*-', 'LineWidth',1)
% hold off
% lg =legend('$g(x)$','$u(x)$');
% lg.FontSize = 16;
% xlabel('$x$')
% title('Convex')
 

% --- Energy

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


%Plots the constraint error
figure(8)
semilogy(iter,error_H_count(1:count-1),'LineWidth',2)
set(gca,'FontSize',16)
hold on
semilogy(count_min, error_H_count(count_min-1),'*','LineWidth',2)
hold off
set(0,'defaulttextInterpreter','latex')
xlabel('Iterations')
ylabel('Log scale')
title('error derivative')




