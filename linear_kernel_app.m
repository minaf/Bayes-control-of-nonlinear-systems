function [c] = linear_kernel_app(med, R, Mz, lx, b)
n = size(med, 1);
z = tf('z');
w = impulse(Mz, n);
w1 = w;
w(1,1) = w(1,1)-1;

G_k = tf(med', [1 zeros(1, n-1)], -1);

Phi_tmp = tf([1 0], [1, b], -1);
phi_tmp = impulse((Mz-1)*Phi_tmp, n);
phi = phi_tmp;
B = impulse((Mz-1)*Phi_tmp*G_k, n);
size(B)
for k=1:lx-1
    Phi_tmp = Phi_tmp/z;
    phi_tmp = impulse((Mz-1)*Phi_tmp, n);
    phi_tmp_new = impulse((Mz-1)*Phi_tmp*G_k, n);
    phi = [phi phi_tmp];
    B = [B phi_tmp_new];
end


A = zeros(lx, lx);

for k=1:lx
     T_phi_til1 = toeplitz(phi(:, k), [phi(1, k) zeros(1, n-1)]);
    for l=1:lx
        T_phi_til2 = toeplitz(phi(:, l), [phi(1, l) zeros(1, n-1)]);
        A(k, l) = trace(T_phi_til1'*T_phi_til2*R);
    end
end

c = -inv(B'*B+A)*B'*w1;
end

