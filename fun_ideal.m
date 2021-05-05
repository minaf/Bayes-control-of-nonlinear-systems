function [suma1] = fun_ideal(x, G, Mz, ts, b, lx)
    K = tf(x', [1 b zeros(1, lx-2)], ts);
    s = impulse(Mz-G*K/(1+G*K), 100);
    suma1 = s'*s;
end
