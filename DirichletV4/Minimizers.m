% Numerically find minimizers for Nonsmooth potentials using Modified 
% Split Bregman 

% Gabriela Jaramillo & Shankar Venkataramani


% E = int (u-u_k)^2/(2h) +W[u'] +V[u] dx  x in D

% Here W[u'] can be chosen to be the convex envelope of
% W1(d) = (d^2-1)^2  "double"
% W2(d) = (d^2-1)^2 if d >= 0 and infty if d<0     "double-half"
% W3(d) = (d^2-1)^2( (d-1)^2-1)^2   "triple"
% W4(d) = random function 

% Also V[u] can be chosen to be
% V1[u] = (u-g(x))^2  "convex potential"
% V2[u] = (u^2 -g(x))^2 "non-convex potential"
% the function g(x) can be chosen in section "Matrices for Gauss Seidel"

% We assume homogeneous Dirichlet BC

% We compute the convex envelop of W(d) using the 'Beneath and Beyond'
% algorithm of Y. Lucet.

% We will assume that the function lives on the nodes, the derivative
% lives on the intervals between the nodes.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc;
global lmb a h

lmb = 0.01; %For constraint (appears in Gauss Seidel iteration)
h =0.01;

potential = 'convex';  % can take values 'non-convex', 'convex'
example = 'double-half'; % options are: 'double', 'double-half', 'triple'

%if 'random' example is selected then:
from_file = 'yes'; %options are: 'yes' if random functions is from file, 'no' if random function is defined here

% If picking potentia= non-convex, you have to pick value of g in section
% Matrices for Gauss-Seidel.

nmx= 2^8;                       % number of nodes
alpha = 0; beta =1;             % end points
dx = (beta-alpha)/(nmx-1);      % grid spacing
xx = (alpha:dx:beta)'; 

u0 = 0.2*ones(size(xx));    %initial guess
u_L=0; u_R=1/2;             %boundary conditions


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Convex Hull using Beneath and Beyond Algorithm
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
    well_random= [d, well];

else
    dd = CX';
    vals = CY';
    well_random = 0;
end


%-- Plot of convex hull and potential well

  figure(1)
  plot(dd,vals,'-o','LineWidth',1)
  set( gca,'FontSize',16)
  hold on
  plot(d,well,'LineWidth',1)
  hold off
  xg =xlabel({'$d$'},'Interpreter','latex');
  xg.FontSize =16;
  lg =legend({'$\overline{W}(d)$','$W(d)$'},'Interpreter','latex');
  lg.FontSize = 16;
 axis([dd(1),dd(end),-0.05,3])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Matrices for Guass Seidel with Dirichlet BC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Possible values of g(x)
%-------------------------------------------
% g = sin(2*pi*xx)/4;
 g = zeros(size(xx));
% g = -1/2*xx;
% g = sin(2*pi*xx)/6+exp(xx)/2;
% g = exp(xx);
% g = sin(4*pi*xx)/12+exp(xx)/2;
% g = -(3/128*16)*(xx).^5 - (xx).^3/12;
% g = -(128/3)*(abs(xx)-0.5).^5 - (abs(xx)-0.5).^3/3;
% g = 1.5*xx;

a = 4*abs(min(g))+0.1;  % for convex splitting (u^2 -g)^2 = 2au^2 +(u^4 -2(g+a)u^2 +g^2)

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

[u,count,Energy, U_min, E_min, count_min, error_H_count] = Split_Bregman_Combined(parameters, vals, dd, g, Lwr, Upr,u0, example,potential);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %   Plots and Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Solution%%%%


%Plot together limiting solution, u,and minimizer, U_min
figure(10)
plot(xx,u,'LineWidth',2)
hold on
plot(xx,U_min,'--','LineWidth',2)
hold off
legend('u last','Umin')

%Calculating derivative
ux = ( u(2:end) - u(1:end-1) )./dx;
Ux = ( U_min(2:end) - U_min(1:end-1)  )./dx;

%Plots derivative of limiting solution, u, and minimizer, U_min
figure(20)
plot(xx(2:end), ux,'LineWidth',2)
hold on
plot(xx(2:end), Ux, 'LineWidth',2)
hold off
legend('u last','U_min')

% Plot derivative of U_min and U_x

% curvature
curve= ( Ux(2:end) - Ux(1:end-1))./dx;
curve_ind = find((curve>10)| (curve<-10));
x_tick = xx(curve_ind);

figure(4)
clf
[hAx,hLine1,hLine2] = plotyy(xx,U_min,xx(2:end),Ux);
legend({'$u(x)$','$u_x(x)$'},'Interpreter','latex')
hLine2.LineStyle = '--';
hLine1.LineWidth = 2;
hLine2.LineWidth = 2;

hAx(1).FontSize = 16;
hAx(2).FontSize = 16;
hAx(1).YColor='k';
hAx(2).YColor ='k';
hAx(2).YGrid ='on';
hAx(2).XGrid ='on';

hAx(2).YTick =[-2 -1 0 1 2];

xlabel('$x$','Interpreter','latex')

ylabel(hAx(1),'$u$','Interpreter','latex') % left y-axis 
ylabel(hAx(2),'$u_x$','Interpreter','latex') % right y-axis


%%% Energy %%%%

M = count-2;
iter = 1:1:count-1;

% %Plot energy at each iteration of Gradient Flow and specifies the
% %iteration with the lowest energy
 
start = count_min-30;
last = count_min+29;

figure(14)
clf
semilogy(iter(1:count-1),Energy(1:count-1),'LineWidth',2)
%plot(iter(1:count-1),Energy(1:count-1),'LineWidth',2)
axis([1 count-1 0.4 1.1])
set(gca,'FontSize',16)
set(0,'defaulttextInterpreter','latex')
hold on 
semilogy(count_min,E_min,'*','LineWidth',2)
xlabel('Iterations')
ylabel('log scale')
title('Energy $\bar{I}$')
axes('Position',[0.4 0.4 .4 .4] )
     box on
     semilogy(iter(start:last),Energy(start:last),'LineWidth',2)
     set(gca,'FontSize',14)
     hold on
     semilogy(count_min-1,E_min,'*','LineWidth',2)
     semilogy(iter(start:last),ones((last-start)+1,1)*Energy(count-1),'--','LineWidth',2)
     hold off
     axis([ start last+2 E_min-10^-6 E_min+10^-6])
 
 
     
%%% Error in derivative %%%%
     
     
 % Plots L^2 for error of derivative (squred)
figure(21)
clf
semilogy(iter,error_H_count(1:count-1),'LineWidth',2)
set(gca,'FontSize',16)
set(0,'defaulttextInterpreter','latex')
hold on
semilogy(count_min,error_H_count(count_min-1),'*','LineWidth',2)
xlabel('Iterations')
ylabel('Log scale')
legend({'$\|D_n -U_n\|_2^2$','$\|D_n -U_{min}\|_2^2$'},'Interpreter','latex')
