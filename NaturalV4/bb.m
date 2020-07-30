function [CX,CY] = bb(X,Y,sm,sp)
%BB implements the 'Beneath and Beyond' algorithm for computing the convex
%envelope. Code is based on work of Y. Lucet.
%
% Input: X,Y with Y(i) = f(X(i)). The X(i) have to be ordered.
% Input: sm,sp are the slopes at the end points. We need sm < 0, sp>0. 
%
% Output: Convex envelope as a piecewise linear function prescribed by
% f**(CVX(j)) = CVY(j)

[~,im] = min(Y-sm*X);
[~,ip] = min(Y-sp*X);

X = X(im:ip);
Y = Y(im:ip);

n = ip-im+1;

if (n <= 2)
    CX=X; CY=Y;
    v=n;
else
    CX(1:2)=[X(1);X(2)];
    CY(1:2)=[Y(1);Y(2)];
	v=2;
    
    for i=3:n
        
        val=((CY(v)-CY(v-1))/(CX(v)-CX(v-1)))*(X(i)-CX(v))+CY(v);
        
        while ((v>2) && (Y(i)<=val)) %Erase points which are not vertices of the convex hull
            v=v-1;
            val=((CY(v)-CY(v-1))/(CX(v)-CX(v-1)))*(X(i)-CX(v))+CY(v);
        end
        
        if v>2
            CX(v+1)=X(i);
            CY(v+1)=Y(i);
            v=v+1;
        else
            if Y(i)> val % Trivial convex hull
                CX(3)=X(i);CY(3)=Y(i);v=3;
            else
                CX(2)=X(i);CY(2)=Y(i);v=2;
            end
        end
    end
end

CX = CX(1:v);
CY = CY(1:v);

    
end

