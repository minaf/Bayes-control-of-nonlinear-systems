function [P, m, c, sigma] = hammerstein_kernel_ident(F, u, y, n, tolerance)

m1 = size(F, 2);
N = size(u, 1);


beta = rand(1);

%%
c = rand(m1, 1);


v = F*c;
sigma = rand(1);

W = toeplitz(v, [v(1) zeros(1, n-1)]);
K = create_TC_kernel(beta, n);

l = 0;
tmp = [];
for t = 1:m1
    tmp_new = toeplitz(F(:, t), [F(1, t), zeros(1, n-1)]);
    tmp = [tmp tmp_new(:)];
end
 tmp = tmp';
while 1

P1 = (K*(W'*W)/sigma^2 + eye(n))\K;
P = 0.5*(P1+P1');
C = P*W'/sigma^2;
m =  C*y;

%ml = (W'*W)\W'*y;

   
    tmp_imp = zeros(m1, m1);
    sve = P+m*m';
    for k = 1:n
        tmp_imp_1 = 0;
        for l = 1:n
            tmp_imp_1 =  tmp_imp_1 + tmp(:, (l-1)*N+1:l*N)*sve(l, k);
        end
        tmp_imp = tmp_imp + tmp_imp_1*tmp(:, (k-1)*N+1:k*N)';
    end
    
    A = tmp_imp;
%     size(F')
%     size(toeplitz([m; zeros(N-n, 1)], [m(1) zeros(1, N-1)])')
%     size(y)
    b = F'*toeplitz([m; zeros(N-n, 1)], [m(1) zeros(1, N-1)])'*y;
    
    
    c_old = c;
    
    c = A\b;
    
    %update sigma 

    v = F*c;
    W = toeplitz(v, [v(1) zeros(1, n-1)]);
    sigma_old = sigma;
    sigma= sqrt(1/N*(norm(y-W*m, 2)^2+trace(W*P*W')));
    
    
    %%
    %update beta
    best_cost = inf;
    beta_old = beta;
    for it = 0.1:0.1:0.9
        cost = (log(det(create_TC_kernel(it, n)))+trace((create_TC_kernel(it, n)))\(P+m*m'));
        %cost = (log(det(create_TC_kernel(it, n)))+ml'/(create_TC_kernel(it, n))*ml);
        if cost<=best_cost
            best_cost = cost;
            beta = it;
        end
    end

%     [lambda, rho] = tools.Q(P+m*m','fir',[n, 0]);
%     lambda = lambda(1);
%     beta = rho(1);

% 

   if norm([sigma; c; beta]-[sigma_old; c_old; beta_old], 2)<=tolerance
       break;
   end
    K = create_TC_kernel(beta, n);
    l = l+1;

end

% z = [y v];
% orders = [0 n 0];
% option = arxRegulOptions('RegularizationKernel','TC');
% arxOpt = arxOptions;
% 
% [Lambda,Reg] = arxRegul(z,orders,option);
% arxOpt.Regularization.Lambda = Lambda;
% arxOpt.Regularization.R = Reg;
% model = arx(z,orders,arxOpt);
% med = getpvec(model, option);
% 
% 
% Phi = toeplitz(v(n:end), v(n:-1:1));
% y_all = y(n:end);
% N = size(y_all, 1);
% 
% % 
% lambda_est = (y_all-Phi*med)'*(y_all-Phi*med)/N;
% K = (eye(n)/(Lambda*Reg))*lambda_est;
% 
% 
% R = K - K*Phi'/(Phi*K*Phi'+lambda_est*eye(size(Phi, 1)))*Phi*K;
% 
% P = (R'+R)/2;
% norm(m-med)
% m = med;

end

