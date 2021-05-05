function [cost] = calculate_cost(x, Mz, G, b, ts)
    lx = size(x, 1);
    K = tf(x',[1 b zeros(1, lx-2)],ts);
    cost = norm(Mz-minreal(K*G/(1+K*G)), 2)^2;
end