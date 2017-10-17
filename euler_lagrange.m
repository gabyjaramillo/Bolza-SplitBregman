function [dydt] = euler_lagrange(~,y)

u = y(1);
p = y(2);

dydt = [p;u./(6*p.^2-2)];

end

