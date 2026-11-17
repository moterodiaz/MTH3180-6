function test_string()
for n = [1 2 3 7 20]
    p = struct('n',n,'M',1.3,'Tf',2,'L',0.9,'c',0.007);
    p.dx = p.L/(n+1);
    p.Uf_func = @(t) 0.002*sin(3*t)+0.001*t^2;
    p.dUfdt_func = @(t) 0.006*cos(3*t)+0.002*t;
    U = (1:n)'*0.001;
    W = -(1:n)'*0.002;
    t = 0.37;
    actual = string_rate_func01(t,[U;W],p);
    expected = zeros(n,1);
    for i = 1:n
        left = 0; right = p.Uf_func(t);
        vl = 0; vr = p.dUfdt_func(t);
        if i > 1
            left = U(i-1); vl = W(i-1);
        end
        if i < n
            right = U(i+1); vr = W(i+1);
        end
        expected(i) = (p.Tf*(left-2*U(i)+right)+p.c*(vl-2*W(i)+vr))/(p.dx*p.M/n);
    end
    assert(norm(actual-[W;expected],inf) < 1e-10);
    [ts,Vs] = string_simulation_template01(p,U,W,[0 0.5],1/2000);
    assert(all(isfinite(Vs),'all') && ts(end) == 0.5);
end
fprintf('Boundary equations and arbitrary forcing: passed (n=1,2,3,7,20)\n');

p = struct('n',3,'M',1,'Tf',1,'L',1,'c',0);
p.Uf_func = @(t) 0*t;
p.dUfdt_func = @(t) 0*t;
errors = zeros(1,3);
for j = 1:3
    h = 0.04/2^(j-1);
    shape = sin((1:3)'*3*pi/4);
    omega = 2*sqrt(12)*sin(3*pi/8);
    [t,V] = string_simulation_template01(p,shape,zeros(3,1),[0 4],h);
    exact = cos(omega*t)*shape';
    errors(j) = max(abs(V(:,1:3)-exact),[],'all');
end
assert(all(errors(1:2)./errors(2:3) > 12));
fprintf('RK4 undamped analytical mode errors: %.3g, %.3g, %.3g\n',errors);

a = string_modes(1/600);
b = string_modes(1/1200);
for k = 1:numel(a)
    s = a(k);
    exact = real(exp(1i*s.omega*s.t)*s.q.');
    err = max(abs(s.V(:,1:s.p.n)-exact),[],'all');
    fine = b(k).V(1:2:end,:);
    change = max(abs(s.V(:,1:s.p.n)-fine(:,1:s.p.n)),[],'all');
    u = [zeros(size(s.t)) s.V(:,1:s.p.n) s.p.Uf_func(s.t)];
    slope = max(abs(diff(u,1,2)/s.p.dx),[],'all');
    eigenfreqs = sqrt(sort(eig((s.p.Tf/(s.p.dx*(s.p.M/s.p.n)))* ...
        (2*eye(s.p.n)-diag(ones(s.p.n-1,1),1)-diag(ones(s.p.n-1,1),-1)))));
    assert(abs(eigenfreqs(s.mode)-s.omega) < 1e-12);
    assert(err < 1e-8 && change < 1e-8 && slope < 0.15);
    fprintf('n=%d mode=%d: omega=%.6f, exact error=%.3g m, half-step change=%.3g m, max slope=%.4f\n', ...
        s.p.n,s.mode,s.omega,err,change,slope);
end
fprintf('All tests passed.\n');
end
