function [controller] = fun_init(trials, Mz, M, ts, b, lx, med, R)
tmp = mvnpdf(med, med, R)/mvnpdf(med, med, 2*R);

z = tf('z');

s = impulse(Mz, M);
res = 0;
out = 0;
for k = 1:size(trials, 2)
    G = tf(trials(:, k)', [1, zeros(1, size(trials, 1)-1)], ts);
    Phi = [];
    for l = 1:lx
        tmp1 = impulse((1-Mz)*G*(1/(1+b*z^(-1)))*z^(-l+1), M);
        Phi = [Phi tmp1];
    end
    tmp1 = mvnpdf(trials(:, k), med, R)/mvnpdf(trials(:, k), med, 2*R)/tmp;
    res = res+(Phi'*Phi);
    out = out+(Phi'*s);
end
    controller = res\out;
end
