function cases = string_modes(h)
if nargin == 0
    h = 1/600;
end
cases = struct([]);
for n = [2 3]
    p = struct('n',n,'M',1,'Tf',1,'L',1,'c',0.004);
    p.dx = p.L/(n+1);
    D = 2*eye(n)-diag(ones(n-1,1),1)-diag(ones(n-1,1),-1);
    for j = 1:n
        omega = 2*sqrt(p.Tf/(p.dx*(p.M/n)))*sin(j*pi/(2*(n+1)));
        Z = (p.Tf+1i*omega*p.c)/p.dx;
        f = zeros(n,1);
        f(end) = Z;
        q = (Z*D-omega^2*(p.M/n)*eye(n))\f;
        A = 0.02/max(abs(q));
        q = A*q;
        p.Uf_func = @(t) A*cos(omega*t);
        p.dUfdt_func = @(t) -A*omega*sin(omega*t);
        % steady-state start
        [t,V] = string_simulation_template01(p,real(q),real(1i*omega*q),[0 8],h);
        k = numel(cases)+1;
        cases(k).t = t;
        cases(k).V = V;
        cases(k).p = p;
        cases(k).mode = j;
        cases(k).omega = omega;
        cases(k).A = A;
        cases(k).q = q;
    end
end
end
