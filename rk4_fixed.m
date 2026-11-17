function [t,V] = rk4_fixed(rate,tspan,V0,h)
assert(h > 0 && isfinite(h),'Step size must be positive and finite.');
assert(numel(tspan) == 2 && tspan(2) > tspan(1),'Use an increasing time interval.');
steps = ceil((tspan(2)-tspan(1))/h);
h = (tspan(2)-tspan(1))/steps;
t = linspace(tspan(1),tspan(2),steps+1)';
V = zeros(steps+1,numel(V0));
V(1,:) = V0(:)';
for i = 1:steps
    y = V(i,:)';
    k1 = rate(t(i),y);
    k2 = rate(t(i)+h/2,y+h*k1/2);
    k3 = rate(t(i)+h/2,y+h*k2/2);
    k4 = rate(t(i)+h,y+h*k3);
    V(i+1,:) = (y+h*(k1+2*k2+2*k3+k4)/6)';
end
end
