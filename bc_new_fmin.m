function [x] = bc_new_fmin(R, med, Mz, lx, b, ts, M, MC)
% z = [y1 u1];
% opt = impulseestOptions('RegularizationKernel','SS');
% sys = impulseest(z, n, opt);
% med = getpvec(sys, opt);
% R = getcov(sys);
% norm(R(1:end-1, 1:end-1))
%n = size(med, 1);
n = size(med, 1);
trials = zeros(n, MC);
R = (R+R')/2;
for k = 1:MC-1
    trials(:, k) = chol(2*R)*randn(size(med, 1),1)+med;
end

trials(:, MC) = med;

p = impulse(Mz);

fu = @(x) fun_bc(x, trials, Mz, M, ts, b, lx, med, R);
c = linear_kernel_app(med, R, Mz, lx, b);
[x, fval]= fmincon(fu, c, [], []);
end

