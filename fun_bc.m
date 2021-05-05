function [suma1] = fun_bc(x, trials, Mz, M, ts, b, lx, med, R)
suma1 = 0;
tmp = mvnpdf(med, med, R)/mvnpdf(med, med, 2*R);
l = (M+1):-1:1;
for k = 1:size(trials, 2)
    K = tf(x',[1 b zeros(1, lx-2)],ts);
    G = tf(trials(:, k)', [1, zeros(1, size(trials, 1)-1)], -1);
    s = impulse(Mz-K*G/(1+K*G), M);
    tmp1 = mvnpdf(trials(:, k), med, R)/mvnpdf(trials(:, k), med, 2*R)/tmp;
    s1 = s.*l';
    suma1 = suma1 + tmp1*(s1'*s);
end
end
