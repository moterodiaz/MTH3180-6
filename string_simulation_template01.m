function [t,V,p] = string_simulation_template01(p,U0,dUdt0,tspan,h)
if nargin == 0
    p = struct('n',3,'M',1,'Tf',1,'L',1,'c',0.004);
    p.Uf_func = @(t) 0.001*cos(2.65*t);
    p.dUfdt_func = @(t) -0.00265*sin(2.65*t);
    U0 = zeros(p.n,1);
    dUdt0 = zeros(p.n,1);
    tspan = [0 8];
    h = 1/600;
end
assert(p.n >= 1 && p.n == fix(p.n),'n must be a positive integer.');
assert(p.M > 0 && p.L > 0 && p.Tf > 0 && p.c >= 0,'Invalid string parameters.');
assert(numel(U0) == p.n && numel(dUdt0) == p.n,'Initial conditions must match n.');
p.dx = p.L/(p.n+1);
[t,V] = rk4_fixed(@(t,y) string_rate_func01(t,y,p),tspan,[U0(:);dUdt0(:)],h);
if nargout == 0
    animate_string(t,V,p,'string_animation.mp4');
end
end
