function [K_ker] = create_TC_kernel(beta, n)
if (beta<0)
    beta = 0;
end

if(beta>1)
    beta = 1;
end
K = zeros(n);
for k = 1:n
    for l = k:n
        K(k, l) = beta^l;
    end
end

K_ker = K + K' - diag(diag(K)); 

end

