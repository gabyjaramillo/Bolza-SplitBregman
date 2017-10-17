%% Global problem. Use analytic and relaxation by method of moments
% This code mixes Split Bregman and Gradient Flow
% E = int (u-u_k)^2/tau^2 + u^2 + (u'^2-1)^2 dx  x \in [0,1] u'>0

% We will assume that the function lives on the nodes, the derivative and
% the Youngs measure lives on the intervals between the nodes.

%function SplitBregmanGF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
nstep = 10000; % maximum number of iterations
nGS = 10; %number of Gauss-Seidel iterations per update of u,v and b.
error = zeros(1,nstep);

err = 1.0;
tol = 1.0e-12;

nmx= 100; % number of nodes
dx = 1/(nmx-1); %grid spacing
xx = (0:dx:1)'; %grid of x values. 

lmb = 2;
tau = 0.001;

%Set up the matrices for Guass Seidel with Dirichlet BC
e = ones(nmx-2,1);
Upr = spdiags(-e,1,nmx-2,nmx-2);
Lwr = spdiags([-e 2*e],-1:0,nmx-2,nmx-2); 

Upr = lmb*Upr/(2*dx^2);
Lwr = (1 + 1/(2* tau))*speye(nmx-2)+lmb*Lwr/(2*dx^2);

% weights to do trapezoid rule for computing the energy
msk = ones(nmx,1); 
msk(1) = 1/sqrt(2.);
msk(end) = 1/sqrt(2.);

% In case we want to include a forcing term
%g0 = sin(2*pi*xx)/4;
g0 = zeros(size(xx));
%g0 = sin(2*pi*xx)/6+exp(xx)/2;
%g0 = exp(xx);
%g0 = sin(4*pi*xx)/12+exp(xx)/2;
%g0 = -(3/128*16)*(xx).^5 - (xx).^3/12;
%g0 = -(128/3)*(abs(xx)-0.5).^5 - (abs(xx)-0.5).^3/3;
%g0 = 1.5*xx;


% If we want to reset u. Otherwise it will use the result of the previous
% computation. Useful if we want to anneal.
reset_flag = true;
if (reset_flag)
  u = g0;
%     u = 0*xx;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Shrink Operator
%  We compute the subgradients at the "breaks" to identify the inverse of
%  the shrink operator

N =60;
deltad =2/(N+1);
daux = (0:deltad:2)';
uu = Obstacle2(N,daux,deltad);
xistar = sqrt(2/3);
breaks = [0 linspace(xistar,3/2,40)];
vals =pchip(daux,uu,breaks);


%vals = ppval(wrelax,breaks);

% In general, breaks and vals are obtained by solving an obstacle problem
slopes = (vals(2:end)-vals(1:end-1))./(breaks(2:end)-breaks(1:end-1));

sub = [-1000 slopes; slopes 1000];
sub = sub(:)';

corners = [breaks;breaks];
corners = corners(:)';

inv_fn = corners + sub/lmb;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initializing variables
b = zeros(nmx-1,1);
v = b;
%vtmp = v;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iteration

count = 1; % number of iterations
objctv =1;
while((err > tol) && (count <= nstep))
    
    uold = u;
    %vold = v;
    %bold = b;
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gauss-Siedel update for u   
    wext = v-b;
    u(1) = 0;
    u(end) = 0.5;
    
    for jj=1:nGS
        frcng = g0(2:end-1) + lmb*(wext(1:nmx-2)-wext(2:nmx-1))/(2*dx)...
            -Upr*u(2:end-1)+lmb/(2*dx^2)*[u(1);zeros(nmx-4,1);u(end)]...
            + uold(2:end-1)/(2*tau);
        u(2:end-1) = Lwr\frcng;
    end

 
    % Update for v
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Using the piecewise linear shrink operator from above 
  
    vtmp = (u(2:end,1)-u(1:end-1,1) )/dx +b;
    v = interp1(inv_fn,corners,vtmp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Using relaxation by hand
%     vtmp = (u(2:end,1)-u(1:end-1,1) )/dx +b;
%      vhelp = sqrt(32/27)/lmb + vtmp;
%      v = max(0,vhelp);
%      ind = find( vhelp>sqrt(2/3)); 
%     
%     %%%  Newton iteration for v for curved wall (analytic relaxation, or by
%     %%% hand)
%     vaux = v;
%      for jj=1:nGS
%          
%         vdelta = vtmp - vaux.* (1+ (4/lmb)*(vaux.^2-1));
%         der = 1 + (4/lmb)*(vaux.^2-1) + (8/lmb)*vaux.^2;
% 
%         vhelp2 = vaux + vdelta./der;
%         vaux = vhelp2;
%      end
% 
%     v(ind) = vaux(ind); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update for b
    b = b-v+(u(2:nmx)-u(1:nmx-1))/dx;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking error   
    
    err = norm(msk.*(u-uold))*sqrt(dx);
    error(count) = err;
    count = count+1;
    
end
% End of iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Analytical solution

options = odeset('Events',@(x,y) endpoint(x,y,u(end)));
y0 = [0;sqrt(2/3)];

sol = ode45(@(x,y) euler_lagrange(x,y),[0,1],y0,options);

x = sol.x;
y = sol.y;

xc = sol.xe;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize data

figure(13)
clf
plot(xx,u,'*');
hold on
plot([0 x+xx(end)-xc],[0 y(1,:)],'-o');
hold off



figure(23)
clf
iters = 1:min(count,nstep);
semilogy(iters,error(iters));
title('Error')


