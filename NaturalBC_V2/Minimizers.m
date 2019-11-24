% Numerically find minimizers for Nonsmooth potentials using Modified 
% Split Bregman 

% Gabriela Jaramillo & Shankar Venkataramani


% E = int (u-u_k)^2/(2h) +W[u'] +V[u] dx  x in D

% Here W[u'] can be chosen to be the convex envelope of
% W1(d) = (d^2-1)^2  "double"
% W2(d) = (d^2-1)^2 if d >= 0 and infty if d<0     "double-half"
% W3(d) = (d^2-1)^2( (d-1)^2-1)^2   "triple"

% Also V[u] can be chosen to be
% V1[u] = (u- g(x))^2   
% V2[u] = (u^2 -g(x))^2  The function g(x) has to be >0 in this case.

% We assume natural BC ( u_x = b+d) at end points

% We compute the convex envelop of W(d) as an obstacle problem

% We will assume that the function lives on the nodes, the derivative
% lives on the intervals between the nodes.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc;
global lmb a h

%a = 4;   %For convex splitting, a>2
lmb = 0.01; %For constraint (appears in Gauss Seidel iteration)
h =0.01;


example = 'double'; % options are: 'double', 'double-half', 'triple'
potential = 'non-convex';  % can take values 'non-convex', 'convex'
% If picking potentia= non-convex, you have to pick value of g in section
% Matrices for Gauss-Seidel. Make sure g>0 on the interval if choosing
% nonconvex

nmx= 2^7;                       % number of nodes
alpha = -1; beta =1;             % end points
dx = (beta-alpha)/(nmx-1);      % grid spacing
xx = (alpha:dx:beta)'; 

u0 = 0.1*ones(size(xx));    % initial guess and BC


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Shrink Operator
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
% g = sin(2*pi*xx)/4 + 1/2;
 g = sin(2*pi*xx)/4;
% g = ones(size(xx));
% g = sin(2*pi*xx)/6+exp(xx)/2;
% g = exp(xx);
% g = sin(4*pi*xx)/12+exp(xx)/2;
% g = -(128/3)*(abs(xx)-0.5).^5 - (abs(xx)-0.5).^3/3;
% g = 1.5*xx;

a = 4*abs(min(g))+0.1;
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

[u,error, count] = Split_Bregman_Combined2(parameters, vals, dd, g, Lwr, Upr,u0);


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ux = ( u(2:end)-u(1:end-1) )./ dx;
        
switch example
    case 'double'
        E1 = max((ux.^2 - 1).^2,0);
        ind = find( abs(ux)<=1);
        E1(ind) = 0;
        
    case 'double-half'
        E1 = max((ux.^2 - 1).^2,0);
        ind = find( abs(ux)<=1);
        E1(ind) = 0;
        
    case 'triple'
        E1 = max((ux.^2 - 1).^2.*((ux-2).^2-1).^2,0);
        ind = find( -1<= ux <=3);
        E1(ind) = 0;
        
end

switch potential
    case 'convex'
        E2 = (u-g).^2;
    case 'non-convex'
        E2 = (u.^2-g).^2;
end

E1 = sum((E1(1:end-1) + E1(2:end))*(dx/2));
E2 = sum((E2(1:end-1) + E2(2:end))*(dx/2));

energy = E1 +E2;
fprintf('Energy is %f\n',energy)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %   Saving file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


filename=sprintf('CS_V2_negativeg.txt'); 
    fileID = fopen(filename,'w');
    fprintf(fileID, '%f ',u);
    fclose(fileID);
   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %   Plots and Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%Plot u and ux
figure(1) % new figure
[hAx,hLine1,hLine2] = plotyy(xx,u,xx(1:end-1),ux);
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
hAx(1).XTick =[-1 -0.66 0 1]; % This needs to change for
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
plot(xx(1:end-1), ux, 'LineWidth',2)
set( gca,'FontSize',16)
set(groot, 'defaultLegendInterpreter','latex');
hold on
plot(xx(1:end-1), dd_sqrtg ,'--', 'LineWidth',2)
hold off
lg =legend('$u_x$','$(\sqrt{g})_x$');
lg.FontSize = 16;


%Plot u and sqrt(g)
figure(3)
plot(xx, sqrt(g) , 'LineWidth',2)
set( gca,'FontSize',16)
set(groot, 'defaultLegendInterpreter','latex');
hold on
plot(xx(1:2:end), u(1:2:end),'*', 'LineWidth',2)
hold off
lg =legend('$\sqrt{g(x)}$','u(x)');
lg.FontSize = 16;
xlabel('$x$')


%%% For convex potential
%%Plot of u and derivative of g

% figure(2)
% plot(xx(1:end-1), ux, 'LineWidth',2)
% set( gca,'FontSize',16)
% hold on
% plot(xx(1:end-1), d_g ,'--', 'LineWidth',2)
% hold off
% lg =legend('$u_x$','$g_x$');
% lg.FontSize = 16;

%%Plot ux and derivative of sqrt(g)
% figure(4) % new figure
% [hAx,hLine1,hLine2] = plotyy(xx(1:end-1),ux,xx(1:end-1),dd_sqrtg);
% set(groot, 'defaultLegendInterpreter','latex');
% set(0,'defaultTextInterpreter','latex')
% lg = legend('$u_x$','$(\sqrt{g})_x$');
% hLine2.LineStyle = '--';
% hLine1.LineWidth = 2;
% hLine2.LineWidth = 2;
% hAx(1).FontSize = 16;
% hAx(2).FontSize = 16;
% hAx(1).YColor='k';
% hAx(2).YColor ='k';
% hAx(1).YGrid ='on';
% hAx(1).XGrid ='on';
% %hAx(1).XTick =[-1 -0.66 0 1]; % This needs to change for
% %different graphs
% xlabel('x')
% ylabel(hAx(1),'$u_x$') % left y-axis 
% ylabel(hAx(2),'$(\sqrt{g})_x$') % right y-axis
% lg.FontSize =16;




 