function dVdt = string_rate_func01(t,V,p)
n = p.n;
U = V(1:n);
dUdt = V(n+1:2*n);

% endpoints
u = [0; U; p.Uf_func(t)];
v = [0; dUdt; p.dUfdt_func(t)];
d2Udt2 = (p.Tf*diff(u,2) + p.c*diff(v,2))/(p.dx*(p.M/n));
dVdt = [dUdt; d2Udt2];
end
